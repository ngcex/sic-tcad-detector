"""Lightweight animation fallback: matplotlib frames -> mp4 (ffmpeg) or gif (Pillow).

Consumes the same ``animation_data(config)`` dict as the manim scenes, so both engines
show the identical real run_cv / run_field data; nothing is interpolated between solved
biases (each solved bias is held for ``hold_s``). Scenes: ``StackScene`` (layers appear
top to bottom, depth to scale) and ``DepletionScene`` (stack bar with the depletion
region, E(x) below, bias / W / peak |E| in the title). mp4 is H.264 yuv420p (plays in
browsers and ``st.video``); gif needs no ffmpeg.
"""

from __future__ import annotations

import shutil
import tempfile
import time
from pathlib import Path

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure

from etna.viz import style as S

FRAME_FORMATS = ("mp4", "gif")


def ffmpeg_available() -> "tuple[bool, str]":
    path = shutil.which("ffmpeg")
    if path:
        return True, ""
    return False, "ffmpeg not found on PATH (needed for mp4; gif works without it)"


def _footer(fig, d, extra=""):
    txt = d["label"] + (f"\n{extra}" if extra else "")
    if d.get("provenance"):
        txt += "\n" + d["provenance"]
    fig.text(0.01, 0.01, txt, fontsize=6, color="#555555", ha="left", va="bottom")


def _stack_frames(fig, d):
    layers = d["layers"]
    total = layers[-1]["bottom_um"]
    ax = fig.add_axes([0.12, 0.14, 0.32, 0.76])
    ax.set_ylim(total, 0)
    ax.set_xlim(0, 1)
    ax.set_xticks([])
    ax.grid(False)
    ax.set_ylabel(S.UNITS["depth"][1])
    ax.minorticks_on()
    ax.tick_params(axis="x", which="both", bottom=False, top=False)
    fig.suptitle("Device stack (depth to scale)")
    _footer(fig, d, "lateral width schematic, not to scale (only the area is known)")
    from matplotlib.patches import Rectangle

    yield  # empty axes first
    for L in layers:
        ax.add_patch(
            Rectangle(
                (0, L["top_um"]),
                1,
                L["thickness_um"],
                facecolor=S.LAYER_COLORS.get(L["key"], "#ccc"),
                edgecolor="k",
                lw=0.6,
            )
        )
        fig.text(
            0.47,
            0.90 - 0.76 * (0.5 * (L["top_um"] + L["bottom_um"]) / total),
            f"{L['name']}: {L['thickness_um']:.3g} µm\n{L['doping']}",
            fontsize=7,
            va="center",
        )
        yield


def _depletion_frames(fig, d):
    layers = d["layers"]
    total = layers[-1]["bottom_um"]
    frames = d["frames"]
    ymax = float(d["E_axis_max"]) * 1e-3
    axb = fig.add_axes([0.09, 0.70, 0.88, 0.12])
    axe = fig.add_axes([0.09, 0.17, 0.88, 0.45])
    from matplotlib.patches import Rectangle

    for L in layers:
        axb.add_patch(
            Rectangle(
                (L["top_um"], 0),
                L["thickness_um"],
                1,
                facecolor=S.LAYER_COLORS.get(L["key"], "#ccc"),
                edgecolor="k",
                lw=0.5,
            )
        )
        axb.text(
            L["top_um"] + L["thickness_um"] / 2,
            -0.15,
            L["key"],
            ha="center",
            va="top",
            fontsize=6,
        )
    axb.set_xlim(0, total)
    axb.set_ylim(0, 1)
    axb.axis("off")
    depl = Rectangle((d["junction_um"], 0), 0, 1, facecolor="#ffd92f", alpha=0.6, lw=0)
    axb.add_patch(depl)
    axe.set_xlim(0, total)
    axe.set_ylim(0, ymax)
    axe.set_xlabel(S.UNITS["depth"][1])
    axe.set_ylabel(S.UNITS["field"][1])
    for xx in (d["junction_um"], d["epi_end_um"]):
        axe.axvline(xx, color="#888", lw=0.6, ls="--")
    (line,) = axe.plot([], [], color=S.PALETTE[0], lw=1.6)
    title = fig.suptitle("")
    rt = d.get("reach_through_V")
    extra = f"W: run_cv ({d['cv_solver']}); E(x): run_field ({frames[0]['solver']})"
    if rt is not None:
        extra += f"; reach-through {rt:.2f} V (Debye-length definition)"
    extra += "; E clipped at the axis top at abrupt doping steps (nm spikes)"
    _footer(fig, d, extra)
    for fr in frames:
        w = max(fr["W_um"], 1e-3)
        depl.set_width(w)
        line.set_data(
            np.asarray(fr["x_um"]),
            np.clip(np.asarray(fr["E_V_per_cm"]) * 1e-3, 0, ymax),
        )
        pk = fr["peak_E_V_per_cm"]
        pk_txt = "n/a" if pk is None else f"{pk * 1e-3:.3g} kV/cm"
        title.set_text(
            f"V = {fr['bias_V']:.1f} V    W = {fr['W_um']:.2f} µm    peak |E| = {pk_txt}"
        )
        yield


_SCENES = {
    "StackScene": (_stack_frames, 1.0),
    "DepletionScene": (_depletion_frames, 0.6),
}


def render_frames(
    data: dict,
    scene: str = "DepletionScene",
    fmt: str = "mp4",
    out_dir=None,
    hold_s: "float | None" = None,
    size_px: "tuple[int, int]" = (960, 540),
    dpi: int = 120,
    fps: int = 25,
) -> "tuple[Path, float]":
    """Render ``scene`` from ``data`` to ``fmt`` (mp4 / gif); (path, wall seconds).

    Each distinct frame (one per solved bias / per layer) is drawn once and held for
    ``hold_s`` seconds (the last one three times as long): gif via Pillow with per-frame
    durations, mp4 via the ffmpeg concat demuxer (H.264, yuv420p, ``fps`` output).
    """
    from matplotlib.backends.backend_agg import FigureCanvasAgg
    from PIL import Image

    if scene not in _SCENES:
        raise ValueError(f"scene must be one of {tuple(_SCENES)}")
    if fmt not in FRAME_FORMATS:
        raise ValueError(f"fmt must be one of {FRAME_FORMATS}")
    if fmt == "mp4":
        ok, why = ffmpeg_available()
        if not ok:
            raise RuntimeError(why)
    if size_px[0] % 2 or size_px[1] % 2:
        raise ValueError("size_px must be even (yuv420p)")
    gen, default_hold = _SCENES[scene]
    hold = default_hold if hold_s is None else float(hold_s)
    out = (
        Path(out_dir)
        if out_dir is not None
        else Path(tempfile.mkdtemp(prefix="etna_frames_"))
    )
    out.mkdir(parents=True, exist_ok=True)
    path = out / f"{scene}.{fmt}"
    t0 = time.perf_counter()
    images = []
    with mpl.rc_context(S.RC):
        fig = Figure(figsize=(size_px[0] / dpi, size_px[1] / dpi), dpi=dpi, layout="none")
        canvas = FigureCanvasAgg(fig)
        for _ in gen(fig, data):
            canvas.draw()
            rgba = Image.frombuffer(
                "RGBA", canvas.get_width_height(), bytes(canvas.buffer_rgba())
            )
            images.append(rgba.convert("RGB"))
    durations = [hold] * len(images)
    durations[-1] = 3 * hold
    if fmt == "gif":
        images[0].save(
            path,
            save_all=True,
            append_images=images[1:],
            duration=[int(1000 * d) for d in durations],
            loop=0,
        )
    else:
        _encode_mp4(images, durations, path, fps)
    if not path.exists() or path.stat().st_size == 0:
        raise RuntimeError(f"frame render of {scene} produced no file")
    return path, time.perf_counter() - t0


def _encode_mp4(images, durations, path: Path, fps: int) -> None:
    import subprocess

    with tempfile.TemporaryDirectory(prefix="etna_frames_png_") as td:
        lines = []
        files = []
        for i, (im, d) in enumerate(zip(images, durations)):
            f = Path(td) / f"f{i:04d}.png"
            im.save(f)
            files.append(f)
            lines += [f"file '{f}'", f"duration {d:.3f}"]
        lines.append(f"file '{files[-1]}'")  # concat demuxer: repeat the last file
        lst = Path(td) / "list.txt"
        lst.write_text("\n".join(lines) + "\n")
        cmd = [
            "ffmpeg",
            "-y",
            "-loglevel",
            "error",
            "-f",
            "concat",
            "-safe",
            "0",
            "-i",
            str(lst),
            "-vf",
            f"fps={int(fps)}",
            "-c:v",
            "libx264",
            "-pix_fmt",
            "yuv420p",
            "-movflags",
            "+faststart",
            str(path),
        ]
        proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)
        if proc.returncode != 0:
            raise RuntimeError(f"ffmpeg failed: {proc.stderr[-1500:]}")

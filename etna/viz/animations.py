"""Two manim scenes rendered from REAL simulator data.

- ``StackScene``: the device stack, depth to scale, with thickness quotes and doping.
- ``DepletionScene``: depletion growth and the E(x) profile vs reverse bias, with the
  bias, depletion width W and the peak field on screen (the peak excludes the nm Debye
  spike at the abrupt epi/n+ interface, ``etna.api.results.peak_field_1d``).

Data: ``animation_data(config)`` runs ``etna.run_cv`` (W per bias) and ``etna.run_field``
(E(x) per bias) at 0..-60 V and returns a JSON-able dict; nothing is synthesised. The
scenes are rendered by a ``python -m manim`` subprocess (``_scenes.py``) reading that
JSON, so manim never runs inside the caller's process (e.g. Streamlit). Low quality
(``-ql``, 480p15) by default. Text only (no LaTeX).
"""

from __future__ import annotations

import json
import logging
import os
import subprocess
import sys
import tempfile
import time
from pathlib import Path

import numpy as np

from etna.export.cross_section import stack_layers
from etna.export.labels import VALIDATION_LABEL

logger = logging.getLogger(__name__)

SCENES = ("StackScene", "DepletionScene")
DEFAULT_BIASES = tuple(float(v) for v in np.linspace(0.0, -60.0, 13))  # 5 V steps
_SCENES_FILE = Path(__file__).with_name("_scenes.py")
_INSTALL_HINT = (
    "install the optional extra with `uv sync --extra dev --extra viz` (macOS: pycairo needs "
    "`brew install pkgconf cairo pango`)"
)


def manim_available() -> "tuple[bool, str]":
    """(True, '') if the manim package is installed, else (False, reason).

    ``importlib.util.find_spec`` only: no import (``import manim`` costs ~9 s cold, and
    Streamlit pages call this on every rerun). manim itself runs only in the render
    subprocess (on click); a broken cairo / pango build shows up there as
    render_animation's RuntimeError, and ``animate(engine="auto")`` then falls back to
    the matplotlib frames.
    """
    import importlib.util

    try:
        spec = importlib.util.find_spec("manim")
        err = "No module named 'manim'"
    except (ImportError, ValueError) as e:  # ValueError: sys.modules['manim'] is None
        spec, err = None, f"{type(e).__name__}: {e}"
    if spec is None:
        return False, f"manim is not importable ({err}); {_INSTALL_HINT}"
    return True, ""


def animation_data(config, biases=DEFAULT_BIASES, max_points: int = 400) -> dict:
    """Real per-bias W and E(x) for ``config`` (1D only), JSON-able.

    W comes from ``etna.run_cv`` on the same bias grid, E(x) from ``etna.run_field`` at
    each bias (node field, downsampled to at most ``max_points`` nodes for drawing; the
    peak is taken on the full profile).
    """
    import etna
    from etna.api.results import peak_field_1d

    if config.half_width_um is not None:
        raise ValueError("animations need a 1D device (set half_width_um to None)")
    biases = [float(b) for b in biases]
    if len(biases) < 2 or any(b > 0 for b in biases):
        raise ValueError("biases: at least two reverse (<= 0 V) values")
    b_min = min(biases)
    # 0.5 V C-V grid (cheap on the batched engine); W is interpolated to each frame bias
    n = int(round(-b_min / 0.5)) + 1
    cv = etna.run_cv(config, 0.0, b_min, max(n, 2))
    W_um = np.asarray(cv.metadata["depletion_widths"], dtype=float) * 1e4
    cv_x = np.abs(np.asarray(cv.x, dtype=float))
    o = np.argsort(cv_x)
    frames = []
    for v in biases:
        f = etna.run_field(config, bias_V=v)
        x = np.asarray(f.x, dtype=float)
        E = np.asarray(f.y, dtype=float)
        pk = peak_field_1d(x, E, config)
        keep = np.unique(
            np.linspace(0, x.size - 1, min(max_points, x.size)).astype(int)
        )
        frames.append(
            dict(
                bias_V=v,
                W_um=float(np.interp(abs(v), cv_x[o], W_um[o])),
                x_um=x[keep].tolist(),
                E_V_per_cm=np.nan_to_num(E[keep]).tolist(),
                peak_E_V_per_cm=None if pk is None else pk[0],
                peak_depth_um=None if pk is None else pk[1],
                peak_excludes_interface=bool(pk and pk[2]),
                solver=f.metadata.get("solver", ""),
            )
        )
    peaks = [fr["peak_E_V_per_cm"] for fr in frames if fr["peak_E_V_per_cm"]]
    from etna.viz.provenance import footer_text, provenance

    prov = provenance(
        config=config,
        solver=f"W: {cv.metadata.get('solver', '')}, E(x): {frames[0]['solver']}",
    )
    return dict(
        label=VALIDATION_LABEL,
        provenance=footer_text(prov),
        config_hash=prov["config_hash"],
        layers=[
            {
                k: L[k]
                for k in (
                    "key",
                    "name",
                    "top_um",
                    "bottom_um",
                    "thickness_um",
                    "doping",
                )
            }
            for L in stack_layers(config)
        ],
        junction_um=config.substrate_thickness_um,
        epi_end_um=config.substrate_thickness_um + config.epi_thickness_um,
        cv_solver=cv.metadata.get("solver", ""),
        reach_through_V=cv.metadata.get("reach_through_V"),
        E_axis_max=1.15 * max(peaks) if peaks else 1.0,
        frames=frames,
    )


def render_animation(
    data: dict,
    scene: str = "DepletionScene",
    out_dir: "str | os.PathLike | None" = None,
    quality: str = "l",
    timeout_s: float = 300.0,
) -> "tuple[Path, float]":
    """Render ``scene`` from ``data`` to mp4 in a manim subprocess; (path, wall seconds).

    ``quality``: manim's -q flag letter (l = 480p15 default, m, h). Raises RuntimeError
    with manim's stderr tail if the render fails, or if manim is not importable.
    """
    if scene not in SCENES:
        raise ValueError(f"scene must be one of {SCENES}")
    ok, why = manim_available()
    if not ok:
        raise RuntimeError(why)
    out = (
        Path(out_dir)
        if out_dir is not None
        else Path(tempfile.mkdtemp(prefix="etna_viz_"))
    )
    out.mkdir(parents=True, exist_ok=True)
    data_file = out / "animation_data.json"
    data_file.write_text(json.dumps(data))
    cmd = [
        sys.executable,
        "-m",
        "manim",
        "render",
        f"-q{quality}",
        "--disable_caching",
        "--media_dir",
        str(out / "media"),
        "-o",
        scene,
        str(_SCENES_FILE),
        scene,
    ]
    env = {**os.environ, "ETNA_VIZ_DATA": str(data_file)}
    t0 = time.perf_counter()
    proc = subprocess.run(
        cmd, env=env, capture_output=True, text=True, timeout=timeout_s, cwd=out
    )
    wall = time.perf_counter() - t0
    hits = sorted((out / "media").rglob(f"{scene}.mp4"))
    if proc.returncode != 0 or not hits:
        tail = (proc.stderr or proc.stdout or "")[-2000:]
        raise RuntimeError(
            f"manim render of {scene} failed (exit {proc.returncode}):\n{tail}"
        )
    return hits[-1], wall


ENGINES = ("auto", "manim", "matplotlib")


def animate(
    data: dict,
    scene: str = "DepletionScene",
    fmt: str = "mp4",
    engine: str = "auto",
    out_dir: "str | os.PathLike | None" = None,
) -> "tuple[Path, float, str]":
    """Render ``scene`` from ``data``; (path, wall seconds, engine used).

    ``engine="auto"``: manim when it is importable and ``fmt`` is mp4, else the
    matplotlib-frames fallback (``etna.viz.frames.render_frames``: mp4 via ffmpeg, gif
    via Pillow). Both read the same ``animation_data`` dict.
    """
    from etna.viz.frames import render_frames

    if engine not in ENGINES:
        raise ValueError(f"engine must be one of {ENGINES}")
    if engine == "manim" and fmt != "mp4":
        raise ValueError("manim renders mp4 only; use engine='matplotlib' for gif")
    use_manim = engine == "manim" or (
        engine == "auto" and fmt == "mp4" and manim_available()[0]
    )
    if use_manim:
        try:
            path, wall = render_animation(data, scene, out_dir=out_dir)
            return path, wall, "manim"
        except RuntimeError:
            if engine == "manim":
                raise
            # installed but not working (e.g. a broken cairo build): frames instead
            logger.warning(
                "manim render failed; using the matplotlib frames", exc_info=True
            )
    path, wall = render_frames(data, scene, fmt=fmt, out_dir=out_dir)
    return path, wall, "matplotlib"

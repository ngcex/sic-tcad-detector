"""manim scenes for etna.viz (run by ``python -m manim``; data from $ETNA_VIZ_DATA).

Only ``Text`` is used for labels (Pango; no LaTeX dependency). Do not import this module
from etna: it imports manim at top level and is meant for the render subprocess.
"""

from __future__ import annotations

import json
import os

import numpy as np
from manim import (
    DOWN,
    LEFT,
    RIGHT,
    UP,
    Axes,
    Create,
    DashedLine,
    DoubleArrow,
    FadeIn,
    GrowFromEdge,
    Line,
    Rectangle,
    Scene,
    Text,
    Transform,
    VGroup,
    VMobject,
    Write,
)

# "buffer": the graded n-/n+ transition layer of the default stack (since 2026-10-18), as
# etna.export.cross_section._COLORS
COLORS = {"p+": "#d6604d", "epi": "#92c5de", "buffer": "#4393c3", "n+": "#2166ac"}
FG = "#ffffff"
LABEL_COLOR = "#f4a582"


def _data() -> dict:
    with open(os.environ["ETNA_VIZ_DATA"]) as fh:
        return json.load(fh)


def _footer(d: dict, extra: str = "") -> VGroup:
    lines = [Text(d["label"], font_size=14, color=LABEL_COLOR, slant="ITALIC")]
    if extra:
        lines.append(Text(extra, font_size=12, color="#bbbbbb"))
    if d.get("provenance"):
        lines.append(Text(d["provenance"], font_size=11, color="#bbbbbb"))
    g = VGroup(*lines).arrange(DOWN, buff=0.06, aligned_edge=LEFT)
    return g.to_edge(DOWN, buff=0.12)


class StackScene(Scene):
    """Stack in scale (depth vertical), thickness quotes and doping per layer."""

    def construct(self):
        d = _data()
        layers = d["layers"]
        total = layers[-1]["bottom_um"]
        H, W, x_left, y_top = 5.4, 3.6, -5.2, 2.7
        s = H / total
        title = Text("Device stack (depth to scale)", font_size=30).to_edge(
            UP, buff=0.25
        )
        self.add(_footer(d))  # validation label on screen from the first frame
        self.play(Write(title), run_time=0.8)
        # depth axis
        ax = Line(
            [x_left - 0.35, y_top, 0], [x_left - 0.35, y_top - H, 0], stroke_width=2
        )
        step = 1.0 if total <= 15 else (2.0 if total <= 30 else 5.0)
        ticks = VGroup()
        dv = 0.0
        while dv <= total + 1e-9:
            y = y_top - dv * s
            ticks.add(
                Line([x_left - 0.45, y, 0], [x_left - 0.35, y, 0], stroke_width=2)
            )
            ticks.add(Text(f"{dv:g}", font_size=14).move_to([x_left - 0.7, y, 0]))
            dv += step
        ylab = (
            Text("depth (µm)", font_size=16)
            .rotate(np.pi / 2)
            .move_to([x_left - 1.15, y_top - H / 2, 0])
        )
        self.play(Create(ax), FadeIn(ticks), FadeIn(ylab), run_time=0.6)
        x_arrow = x_left + W + 0.35
        for L in layers:
            y0 = y_top - L["top_um"] * s
            h = L["thickness_um"] * s
            rect = Rectangle(
                width=W,
                height=h,
                fill_color=COLORS[L["key"]],
                fill_opacity=1.0,
                stroke_color=FG,
                stroke_width=1,
            )
            rect.move_to([x_left + W / 2, y0 - h / 2, 0])
            arrow = DoubleArrow(
                [x_arrow, y0, 0],
                [x_arrow, y0 - h, 0],
                buff=0,
                stroke_width=2,
                tip_length=min(0.12, h / 2.5),
                max_tip_length_to_length_ratio=0.45,
            )
            name = Text(
                f"{L['name']}: {L['thickness_um']:g} µm", font_size=22, weight="BOLD"
            )
            dop = Text(L["doping"], font_size=13)
            lab = VGroup(name, dop).arrange(DOWN, buff=0.06, aligned_edge=LEFT)
            lab.next_to([x_arrow + 0.2, y0 - h / 2, 0], RIGHT, buff=0)
            avail = 6.9 - (x_arrow + 0.2)  # frame half-width is ~7.1
            if lab.width > avail:
                lab.scale_to_fit_width(avail)
                lab.next_to([x_arrow + 0.2, y0 - h / 2, 0], RIGHT, buff=0)
            if h < 0.4:  # thin layer (p+): keep its label clear of the next one
                lab.shift(UP * 0.18)
            self.play(GrowFromEdge(rect, UP), Create(arrow), FadeIn(lab), run_time=0.7)
        note = Text(
            "lateral width schematic, not to scale (only the area is known)",
            font_size=14,
            color="#bbbbbb",
        ).next_to([x_left, y_top - H, 0], DOWN, buff=0.15, aligned_edge=LEFT)
        self.play(FadeIn(note), run_time=0.6)
        self.wait(1.5)


class DepletionScene(Scene):
    """Depletion width and E(x) vs reverse bias, from run_cv / run_field."""

    def construct(self):
        d = _data()
        layers = d["layers"]
        total = layers[-1]["bottom_um"]
        frames = d["frames"]
        x0, width = -6.0, 12.0
        sx = width / total
        bar_y, bar_h = 1.75, 0.6
        ymax = float(d["E_axis_max"])

        title = Text("Depletion growth and E(x) vs reverse bias", font_size=28).to_edge(
            UP, buff=0.2
        )
        bar = VGroup()
        for L in layers:
            r = Rectangle(
                width=L["thickness_um"] * sx,
                height=bar_h,
                fill_color=COLORS[L["key"]],
                fill_opacity=1.0,
                stroke_color=FG,
                stroke_width=1,
            )
            r.move_to([x0 + (L["top_um"] + L["thickness_um"] / 2) * sx, bar_y, 0])
            bar.add(r)
        names = VGroup(
            *[
                Text(L["name"], font_size=14).next_to(
                    [
                        x0 + (L["top_um"] + L["thickness_um"] / 2) * sx,
                        bar_y - bar_h / 2,
                        0,
                    ],
                    DOWN,
                    buff=0.08,
                )
                for L in layers
            ]
        )
        names[0].shift(LEFT * 0.25)

        axes = Axes(
            x_range=[0, total, 1],
            y_range=[0, ymax, ymax / 4],
            x_length=width,
            y_length=3.3,
            tips=False,
            axis_config={"include_numbers": False, "stroke_width": 2},
        )
        o = axes.c2p(0, 0)
        axes.shift(np.array([x0 - o[0], -2.65 - o[1], 0]))
        step = 1.0 if total <= 15 else (2.0 if total <= 30 else 5.0)
        xt = VGroup(
            *[
                Text(f"{v:g}", font_size=13).next_to(axes.c2p(v, 0), DOWN, buff=0.08)
                for v in np.arange(0, total + 1e-9, step)
            ]
        )
        yt = VGroup(
            *[
                Text(f"{v:.1e}", font_size=12).next_to(axes.c2p(0, v), LEFT, buff=0.08)
                for v in np.linspace(0, ymax, 5)[1:]
            ]
        )
        xlab = (
            Text("depth (µm)", font_size=15)
            .next_to(axes.c2p(total, 0), DOWN, buff=0.3)
            .shift(LEFT * 0.6)
        )
        ylab = Text("|E| (V/cm)", font_size=15).next_to(
            axes.c2p(0, ymax), UP, buff=0.08
        )
        guides = VGroup(
            *[
                DashedLine(
                    [x0 + xx * sx, bar_y - bar_h / 2, 0],
                    axes.c2p(xx, 0),
                    stroke_width=1,
                    color="#888888",
                )
                for xx in (d["junction_um"], d["epi_end_um"])
            ]
        )
        rt = d.get("reach_through_V")
        extra = f"W: run_cv ({d['cv_solver']}); E(x): run_field ({frames[0]['solver']})"
        if rt is not None:
            extra += f"; reach-through {rt:.2f} V (Debye-length definition)"

        def depl(fr):
            w = max(fr["W_um"], 1e-3)
            r = Rectangle(
                width=w * sx,
                height=bar_h,
                fill_color="#ffd92f",
                fill_opacity=0.55,
                stroke_width=0,
            )
            r.move_to([x0 + (d["junction_um"] + w / 2) * sx, bar_y, 0])
            return r

        def curve(fr):
            x = np.asarray(fr["x_um"])
            e = np.clip(np.asarray(fr["E_V_per_cm"]), 0, ymax)
            c = VMobject(color="#ffd92f", stroke_width=3)
            c.set_points_as_corners([axes.c2p(a, b) for a, b in zip(x, e)])
            return c

        def info(fr):
            pk = fr["peak_E_V_per_cm"]
            pk_txt = "n/a" if pk is None else f"{pk:.3g} V/cm"
            if fr["peak_excludes_interface"]:
                pk_txt += " (excl. spikes at abrupt doping steps)"
            t = Text(
                f"V = {fr['bias_V']:.0f} V    W = {fr['W_um']:.2f} µm    peak |E| = {pk_txt}",
                font_size=19,
            )
            return t.next_to(title, DOWN, buff=0.18)

        clipped = any(max(fr["E_V_per_cm"]) > ymax for fr in frames)
        self.add(_footer(d, extra))  # validation label on screen from the first frame
        self.play(Write(title), FadeIn(bar), FadeIn(names), run_time=0.8)
        self.play(
            Create(axes),
            FadeIn(xt),
            FadeIn(yt),
            FadeIn(xlab),
            FadeIn(ylab),
            FadeIn(guides),
            run_time=0.8,
        )
        if clipped:
            note = Text(
                "E clipped at the axis top at the epi/n+ step (nm Debye spike)",
                font_size=12,
                color="#bbbbbb",
            ).move_to(axes.c2p(0.6 * total, 0.95 * ymax))
            self.add(note)
        cur_d, cur_c, cur_i = depl(frames[0]), curve(frames[0]), info(frames[0])
        self.play(FadeIn(cur_d), Create(cur_c), FadeIn(cur_i), run_time=0.8)
        for fr in frames[1:]:
            self.play(
                Transform(cur_d, depl(fr)),
                Transform(cur_c, curve(fr)),
                Transform(cur_i, info(fr)),
                run_time=0.6,
            )
            self.wait(0.15)
        self.wait(1.5)

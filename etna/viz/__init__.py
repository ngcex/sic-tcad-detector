"""Scientific figures and animations of real etna results.

- Figures (matplotlib, always available): ``figure_for(result, kind)`` -> Figure in the
  publication style with units, validity bands and a provenance footer;
  ``export_figure(fig, "png" | "svg" | "pdf", dpi, size)`` -> bytes; ``kinds_for(obj)``.
- Animations: ``animation_data(config)`` (real run_cv / run_field data, pure numpy +
  etna solvers), then ``render_animation(data, scene)`` (manim, optional extra ``viz``:
  ``uv sync --extra dev --extra viz``; macOS: pycairo needs ``brew install pkgconf cairo
  pango``) or ``render_frames(data, scene, fmt)`` (matplotlib frames -> mp4 via ffmpeg or
  gif via Pillow, no manim needed); ``animate(data, scene, fmt, engine="auto")`` picks one.
- Validation variants: ``etna.viz.validation`` (tilted-alpha shape, W creep).
- 2D carrier maps: ``etna.viz.maps.carrier_state_2d`` (engine2d dark state with n, p).
"""

from etna.viz.animations import (
    SCENES,
    animate,
    animation_data,
    manim_available,
    render_animation,
)
from etna.viz.figures import (
    KINDS,
    export_figure,
    figure_for,
    kinds_for,
    load_let_table,
    save_figure,
)
from etna.viz.frames import ffmpeg_available, render_frames
from etna.viz.provenance import config_hash

__all__ = [
    "SCENES",
    "KINDS",
    "animate",
    "animation_data",
    "config_hash",
    "export_figure",
    "ffmpeg_available",
    "figure_for",
    "kinds_for",
    "load_let_table",
    "manim_available",
    "render_animation",
    "render_frames",
    "save_figure",
]

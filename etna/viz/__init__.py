"""Optional animations of real simulator data (manim, extra ``viz``).

Install with ``uv sync --extra viz`` (macOS: pycairo needs ``brew install pkgconf cairo
pango``). Without manim, ``manim_available()`` returns ``(False, reason)`` and the rest
of etna is unaffected; ``animation_data`` (pure numpy + etna solvers) works either way.
"""

from etna.viz.animations import (
    SCENES,
    animation_data,
    manim_available,
    render_animation,
)

__all__ = ["SCENES", "animation_data", "manim_available", "render_animation"]

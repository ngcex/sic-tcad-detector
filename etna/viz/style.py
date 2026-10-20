"""Publication style for etna.viz figures: one rc set, one palette, fixed sizes, units.

Applied with ``matplotlib.rc_context(RC)`` around figure creation and saving only, never
to the global rcParams (the notebooks and ``etna.core.plotting`` keep their own look).
Text stays editable in vector output (``pdf.fonttype`` 42, ``svg.fonttype`` "none").
"""

from __future__ import annotations

MM_PER_INCH = 25.4

# Figure sizes (width, height) in mm: a journal single column (~89 mm, e.g. NIM A / Phys.
# Med.), a double column (183 mm) and a 16:9 slide.
SIZES_MM = {
    "single": (89.0, 70.0),
    "double": (183.0, 105.0),
    "slide": (254.0, 143.0),
}
DEFAULT_SIZE = "single"
DEFAULT_DPI = 300
FORMATS = ("png", "svg", "pdf")

# Okabe-Ito colour-blind-safe palette (order used for series)
PALETTE = (
    "#0072B2",  # blue
    "#D55E00",  # vermillion
    "#009E73",  # green
    "#CC79A7",  # purple
    "#E69F00",  # orange
    "#56B4E9",  # sky
    "#000000",  # black
)
DATA_COLOR = "#000000"  # measured / digitised points
INVALID_COLOR = "#B0B0B0"  # outside-validity shading
LAYER_COLORS = {  # same keys as etna.export.cross_section.stack_layers
    "p+": "#d6604d",
    "epi": "#92c5de",
    "buffer": "#4393c3",
    "n+": "#2166ac",
}

BASE_FONT_PT = 8.0
FOOTER_FONT_PT = 4.6

RC = {
    "font.family": "sans-serif",
    "font.sans-serif": ["DejaVu Sans", "Arial", "Helvetica"],
    "font.size": BASE_FONT_PT,
    "axes.titlesize": BASE_FONT_PT + 0.5,
    "axes.labelsize": BASE_FONT_PT,
    "axes.linewidth": 0.7,
    "axes.prop_cycle": None,  # filled below (needs cycler)
    "axes.grid": True,
    "grid.linewidth": 0.35,
    "grid.alpha": 0.45,
    "xtick.direction": "in",
    "ytick.direction": "in",
    "xtick.top": True,
    "ytick.right": True,
    "xtick.minor.visible": True,
    "ytick.minor.visible": True,
    "xtick.labelsize": BASE_FONT_PT - 1,
    "ytick.labelsize": BASE_FONT_PT - 1,
    "xtick.major.width": 0.6,
    "ytick.major.width": 0.6,
    "lines.linewidth": 1.2,
    "lines.markersize": 3.0,
    "legend.fontsize": BASE_FONT_PT - 1,
    "legend.frameon": False,
    "pdf.fonttype": 42,
    "ps.fonttype": 42,
    "svg.fonttype": "none",
    "savefig.dpi": DEFAULT_DPI,
    "image.cmap": "viridis",
    "axes.formatter.use_mathtext": True,
    "mathtext.default": "regular",
}


def _finish_rc() -> dict:
    from cycler import cycler

    rc = dict(RC)
    rc["axes.prop_cycle"] = cycler(color=list(PALETTE))
    return rc


RC = _finish_rc()


def size_inches(size) -> "tuple[float, float]":
    """``"single" | "double" | "slide" | (w_mm, h_mm)`` -> (w_in, h_in)."""
    if isinstance(size, str):
        if size not in SIZES_MM:
            raise ValueError(f"size must be one of {sorted(SIZES_MM)} or (w_mm, h_mm)")
        w, h = SIZES_MM[size]
    else:
        try:
            w, h = (float(v) for v in size)
        except (TypeError, ValueError) as e:
            raise ValueError("size must be a preset name or (w_mm, h_mm)") from e
        if not (20.0 <= w <= 600.0 and 20.0 <= h <= 600.0):
            raise ValueError("size in mm must lie in 20..600 for both sides")
    return w / MM_PER_INCH, h / MM_PER_INCH


# Units shown on axes (conversion factor from the SimResult unit, label)
UNITS = {
    "capacitance": (1e12, "C (pF)"),  # F -> pF
    "inv_c2": (1.0, r"$1/C^2$ (F$^{-2}$)"),
    "width_um": (1e4, r"W ($\mu$m)"),  # cm -> um
    "field": (1e-3, r"|E| (kV/cm)"),  # V/cm -> kV/cm
    "potential": (1.0, r"$\psi$ (V)"),
    "current": (1.0, "|I| (A)"),
    "cce": (1.0, "CCE"),
    "depth": (1.0, r"Depth ($\mu$m)"),
    "lateral": (1.0, r"Lateral x ($\mu$m)"),
    "reverse_bias": (-1.0, r"Reverse bias $-V$ (V)"),  # bias V (negative) -> -V
    "temperature": (1.0, "T (K)"),
    "dose_rate": (1.0, "Dose rate (Gy/s)"),
    "fluence": (1.0, r"Proton fluence (cm$^{-2}$)"),
    "let": (1.0, r"LET (keV/$\mu$m)"),
    "lineal": (1.0, r"Lineal energy y (keV/$\mu$m)"),
}

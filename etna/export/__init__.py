"""Engineering exports of a simulated device: spec sheet, cross-section, mesh/profile.

Every function here is pure: it takes a ``DeviceConfig`` and already computed results
(``etna.run_all`` Report, ``etna.run_field`` SimResult) and never runs a simulation.
Every export carries ``VALIDATION_LABEL`` verbatim.

- ``spec_sheet``: ``build_spec`` -> dict; ``spec_to_json`` / ``spec_to_csv`` /
  ``spec_to_markdown`` / ``spec_to_pdf`` (PDF drawn with matplotlib, already a dependency).
- ``cross_section``: dimensioned stack drawing, vertical axis to scale (SVG / PNG).
- ``mesh_export``: doping / potential / field / n / p per node (CSV, legacy-ASCII VTK,
  XML VTU), 1D and 2D.

- ``layout_params``: the user-editable layout / mask dimension file
  (``data/foundry/layout_params.toml``; argument or $ETNA_LAYOUT_PARAMS for another
  one), validated; read by ``gds``, the spec sheet and the cross-section. Dimensions
  still classed PLACEHOLDER are listed in every output and warned about
  (``PlaceholderDimensionWarning``).
- ``gds``: labelled reference layouts ("reference layout, not TCAD-validated"; top cell
  ``*REFERENCE-NOT-FAB`` while any PLACEHOLDER is left), GDSII (optional ``gdstk``) +
  SVG: active p+ square, guard ring(s), corner pad, die outline, dicing street.
"""

from etna.export.labels import VALIDATION_LABEL
from etna.export.spec_sheet import (
    SPEC_SECTIONS,
    build_spec,
    spec_to_csv,
    spec_to_json,
    spec_to_markdown,
    spec_to_pdf,
)
from etna.export.cross_section import cross_section_png, cross_section_svg, stack_layers
from etna.export.mesh_export import profile_to_csv, profile_to_vtk, profile_to_vtu
from etna.export.layout_params import (
    DEFAULT_PARAMS_PATH as LAYOUT_PARAMS_PATH,
    LayoutParams,
    LayoutParamsError,
    PlaceholderDimensionWarning,
    load_layout_params,
    parse_layout_params,
)
from etna.export.gds import (
    DISCLAIMER as LAYOUT_DISCLAIMER,
    gds_available,
    layout_gds,
    layout_svg,
)

__all__ = [
    "VALIDATION_LABEL",
    "SPEC_SECTIONS",
    "build_spec",
    "spec_to_json",
    "spec_to_csv",
    "spec_to_markdown",
    "spec_to_pdf",
    "stack_layers",
    "cross_section_svg",
    "cross_section_png",
    "profile_to_csv",
    "profile_to_vtk",
    "profile_to_vtu",
    "LAYOUT_PARAMS_PATH",
    "LayoutParams",
    "LayoutParamsError",
    "PlaceholderDimensionWarning",
    "load_layout_params",
    "parse_layout_params",
    "LAYOUT_DISCLAIMER",
    "gds_available",
    "layout_gds",
    "layout_svg",
]

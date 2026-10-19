"""Engineering exports of a simulated device: spec sheet, cross-section, mesh/profile.

Every function here is pure: it takes a ``DeviceConfig`` and already computed results
(``etna.run_all`` Report, ``etna.run_field`` SimResult) and never runs a simulation.
Every export carries ``VALIDATION_LABEL`` verbatim.

- ``spec_sheet``: ``build_spec`` -> dict; ``spec_to_json`` / ``spec_to_csv`` /
  ``spec_to_markdown`` / ``spec_to_pdf`` (PDF drawn with matplotlib, already a dependency).
- ``cross_section``: dimensioned stack drawing, vertical axis to scale (SVG / PNG).
- ``mesh_export``: doping / potential / field per node (CSV, legacy-ASCII VTK), 1D and 2D.

- ``gds``: labelled SCHEMATIC reference layouts (top cell ``*REFERENCE-NOT-FAB``), GDSII
  (optional ``gdstk``) + SVG, from the documented / measured-from-figure dimensions only.
  Not a fabrication layout: guard ring, pad and dicing dimensions are not available and
  are not invented.
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
from etna.export.mesh_export import profile_to_csv, profile_to_vtk
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
    "LAYOUT_DISCLAIMER",
    "gds_available",
    "layout_gds",
    "layout_svg",
]

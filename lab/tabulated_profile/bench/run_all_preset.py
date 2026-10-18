"""run_all on the tabulated measured-C-V preset, every section incl. the opt-in 2D sweep,
then the spec sheet (Markdown) of the report.

usage: uv run python lab/tabulated_profile/bench/run_all_preset.py
"""

from __future__ import annotations

import etna
from etna.api.run_all import SECTION_NAMES
from etna.export.spec_sheet import build_spec, spec_to_markdown


def main():
    cfg = etna.measured_cv_tabulated_config(area_cm2=0.04)
    rep = etna.run_all(cfg, sections=SECTION_NAMES)
    for r in rep.rows():
        print("RES", r)
    print(f"RES wall {rep.wall_s:.1f} s, ok {rep.ok}")
    print(spec_to_markdown(build_spec(cfg, rep)))


if __name__ == "__main__":
    main()

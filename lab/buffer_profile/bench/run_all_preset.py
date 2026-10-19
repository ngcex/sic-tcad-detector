"""run_all on the measured-C-V preset, every section incl. the opt-in 2D sweep.

usage: uv run python lab/buffer_profile/bench/run_all_preset.py
"""

from __future__ import annotations

import etna
from etna.api.run_all import SECTION_NAMES


def main():
    cfg = etna.measured_cv_config(area_cm2=0.04)
    rep = etna.run_all(cfg, sections=SECTION_NAMES)
    for r in rep.rows():
        print("RES", r)
    print(f"RES wall {rep.wall_s:.1f} s, ok {rep.ok}")


if __name__ == "__main__":
    main()

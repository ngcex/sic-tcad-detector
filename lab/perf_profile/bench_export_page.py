"""What the Export page computes on EVERY rerun (no button clicked): cross-section
PNG + SVG, gds_available(), layout_gds per device. First call vs repeat.

usage: uv run python lab/perf_profile/bench_export_page.py
"""

from __future__ import annotations

import time

from etna import DeviceConfig
from etna import export as ex

cfg = DeviceConfig()


def timed(label, fn):
    t = time.perf_counter()
    out = fn()
    t1 = time.perf_counter() - t
    t = time.perf_counter()
    fn()
    t2 = time.perf_counter() - t
    print(f"{label:34s} first {t1:7.3f} s  repeat {t2:7.3f} s")
    return out


timed("cross_section_png", lambda: ex.cross_section_png(cfg))
timed("cross_section_svg", lambda: ex.cross_section_svg(cfg))
ok, why = timed("gds_available", ex.gds_available)
print("gds available:", ok, why)
if ok:
    import inspect

    import app.workflows.export as page

    src = inspect.getsource(page._layout_block)
    for dev in ("dosimeter", "microdosimeter"):
        if dev in src:
            timed(f"layout_gds({dev})", lambda d=dev: ex.layout_gds(d))

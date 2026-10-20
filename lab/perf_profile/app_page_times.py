"""Streamlit page load times (AppTest, no button clicks): script-run wall per page.

The first run of the process includes `import etna` and the app imports; later runs
are the per-interaction rerun cost (Streamlit reruns the whole page script on every
widget change). One fresh AppTest per page, run twice.

usage: uv run python lab/perf_profile/app_page_times.py
"""

from __future__ import annotations

import tempfile
import time
from pathlib import Path

from streamlit.testing.v1 import AppTest

PAGES = [
    "home",
    "run_all",
    "cv",
    "cce",
    "field_map",
    "radiation_damage",
    "dark_current",
    "microdosimetry",
    "batch_sweep",
    "export",
]

t = time.perf_counter()
at = AppTest.from_file("app/main.py", default_timeout=600)
at.run()
print(
    f"main.py (sidebar + home), first run in process: {time.perf_counter() - t:.2f} s, exceptions {len(at.exception)}"
)
t = time.perf_counter()
at.run()
print(f"main.py rerun: {time.perf_counter() - t:.2f} s")

SRC = """
import sys
sys.path.insert(0, {root!r})
from app.components.device_sidebar import render_device_sidebar
from app.workflows.{m} import render
render_device_sidebar()
render()
"""
ROOT = str(Path(__file__).resolve().parents[2])
TMP = Path(tempfile.mkdtemp(prefix="etna_app_pages_"))
for m in PAGES:
    script = TMP / f"page_{m}.py"
    script.write_text(SRC.format(m=m, root=ROOT))
    at = AppTest.from_file(str(script), default_timeout=600)
    t = time.perf_counter()
    at.run()
    t1 = time.perf_counter() - t
    t = time.perf_counter()
    at.run()
    t2 = time.perf_counter() - t
    print(f"{m:18s} first {t1:.3f} s, rerun {t2:.3f} s, exceptions {len(at.exception)}")

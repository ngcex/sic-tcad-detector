"""`import etna` stays free of the plotting / statistics stacks (import-time budget).

matplotlib.pyplot (~0.6 s cpu) and scipy.stats (~0.9 s cpu) used to be imported at
module level by etna.core.cv_analysis / microdosimetry / temperature_sweep only for
one plot helper or one linregress call; every process (xdist worker, app start,
notebook kernel, spawn-pool worker) paid them. They are imported inside the
functions now (A/B: lab/perf_profile/ab_import.py). A fresh interpreter is needed: pytest
itself may already have imported them.
"""

from __future__ import annotations

import subprocess
import sys


def test_import_etna_does_not_load_pyplot_or_scipy_stats():
    code = (
        "import sys, etna; "
        "print('LOADED=' + ','.join(m for m in ('matplotlib.pyplot', 'scipy.stats') if m in sys.modules))"
    )
    out = subprocess.run(
        [sys.executable, "-c", code], capture_output=True, text=True, check=True
    )
    # devsim prints its math-library banner on stdout at import: read the marked line
    (line,) = [ln for ln in out.stdout.splitlines() if ln.startswith("LOADED=")]
    loaded = line[len("LOADED=") :]
    assert loaded == "", f"import etna loaded: {loaded}"


def test_lazy_helpers_still_work():
    import numpy as np

    from etna.core.temperature_sweep import extract_temperature_coefficient

    r = extract_temperature_coefficient(
        np.array([250.0, 300.0, 350.0]), np.array([1.0, 2.0, 3.0])
    )
    assert abs(r["slope"] - 0.02) < 1e-15

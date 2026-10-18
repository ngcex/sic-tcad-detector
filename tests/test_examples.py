"""Smoke test: the scripts in examples/ run on the current API (audit-v7 D-10).

Each example runs in a subprocess with cwd = tmp_path and ETNA_FIGURE_DIR = tmp_path,
so it cannot write into the tracked figures/ directory. The example raises on its own
checks (no truncation, monotonic C), so a zero exit status is the assertion.
"""

import os
import pathlib
import subprocess
import sys

import pytest

EXAMPLES = pathlib.Path(__file__).resolve().parents[1] / "examples"


@pytest.mark.parametrize("script", ["cv_example.py"])
def test_example_runs(script, tmp_path):
    env = {**os.environ, "ETNA_FIGURE_DIR": str(tmp_path)}
    proc = subprocess.run(
        [sys.executable, str(EXAMPLES / script)],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        timeout=300,
    )
    assert proc.returncode == 0, proc.stdout[-2000:] + proc.stderr[-2000:]
    assert "truncated=False" in proc.stdout
    assert "decreasing monotonically" in proc.stdout


def test_every_example_is_smoke_tested():
    tested = {"cv_example.py"}
    assert {p.name for p in EXAMPLES.glob("*.py")} == tested

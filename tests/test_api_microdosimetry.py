"""Fast data-pipeline test for the run_microdosimetry facade (LIB-06, bucket-c).

run_microdosimetry is the only bucket-c facade — a pure data pipeline that
loads MC events, aggregates per-event collected energy, and computes an
ICRU-36 lineal energy spectrum. It touches NO devsim, so this test runs the
facade end-to-end against the committed synthetic fixture
(data/synthetic_mc_events.csv) and asserts the SimResult contract.

The CSV path is resolved relative to the repo root via
Path(__file__).parent.parent so the test is CWD-independent (it does not
assume pytest was invoked from the repo root).
"""

from pathlib import Path

import numpy as np

from etna import DeviceConfig, run_microdosimetry
from etna.api.results import SimResult

# CWD-independent path to the committed synthetic MC-events fixture.
MC_CSV_PATH = Path(__file__).parent.parent / "data" / "synthetic_mc_events.csv"


def test_run_microdosimetry_data_pipeline():
    """run_microdosimetry produces a valid y*d(y) SimResult from the fixture."""
    assert MC_CSV_PATH.exists(), f"fixture missing: {MC_CSV_PATH}"

    result = run_microdosimetry(
        DeviceConfig(),
        mc_csv_path=str(MC_CSV_PATH),
        sv_thickness_um=10,
        sv_width_um=150,
    )

    # Output type and sim_type.
    assert isinstance(result, SimResult)
    assert result.sim_type == "microdosimetry"

    # Output shape: x (bin centers) and y (y*d(y)) are equal-length, non-empty.
    x = np.asarray(result.x, dtype=float)
    y = np.asarray(result.y, dtype=float)
    assert len(x) == len(y)
    assert len(x) > 0

    # All values finite (no NaN/inf in the spectrum).
    assert np.all(np.isfinite(x))
    assert np.all(np.isfinite(y))

    # metadata carries the microdosimetric means; Jensen's inequality
    # requires y_D >= y_F for a well-formed lineal energy spectrum.
    assert "y_F" in result.metadata
    assert "y_D" in result.metadata
    assert result.metadata["y_D"] >= result.metadata["y_F"]


def test_sv_width_does_not_enter_l_bar():
    """audit-v7 D-15: with width alone mean_chord_length is the slab value 2 t, so the
    facade default (150 um) and DESIGN-2's 100 um SV give the same l_bar."""
    from etna.core.microdosimetry import mean_chord_length

    default = run_microdosimetry(DeviceConfig(), mc_csv_path=str(MC_CSV_PATH))
    design2 = run_microdosimetry(
        DeviceConfig(), mc_csv_path=str(MC_CSV_PATH), sv_width_um=100.0
    )
    assert default.metadata["l_bar_um"] == design2.metadata["l_bar_um"] == 20.0
    assert mean_chord_length(10.0, sv_width_um=100.0) == 20.0
    np.testing.assert_array_equal(default.y, design2.y)

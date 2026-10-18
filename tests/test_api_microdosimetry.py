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
import pytest

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


def test_in_memory_stream_gives_the_same_spectrum_as_the_path():
    """audit C-17: a file-like object is accepted (the app passes io.BytesIO)."""
    import io

    by_path = run_microdosimetry(DeviceConfig(), mc_csv_path=str(MC_CSV_PATH))
    by_stream = run_microdosimetry(
        DeviceConfig(), mc_csv_path=io.BytesIO(MC_CSV_PATH.read_bytes())
    )
    assert by_stream.metadata["y_F"] == by_path.metadata["y_F"]
    assert by_stream.metadata["y_D"] == by_path.metadata["y_D"]
    np.testing.assert_array_equal(by_stream.y, by_path.y)


@pytest.mark.parametrize(
    "csv, kwargs, match",
    [
        (
            b"event_id,x,y,z,edep\n0,0,0,0,100\n",
            {"sv_thickness_um": 0.0},
            "sv_thickness_um",
        ),
        (
            b"event_id,x,y,z,edep\n0,0,0,0,100\n",
            {"sv_width_um": float("nan")},
            "sv_width_um",
        ),
        (b"event_id,x,y,z,edep\n", {}, "no events"),
        (b"event_id,x,y,z,edep\n0,0,0,0,100\n1,0,0,0,-5\n", {}, "negative"),
        (b"event_id,x,y,z,edep\n0,0,0,0,0\n1,0,0,0,0\n", {}, "deposits any energy"),
    ],
)
def test_meaningless_inputs_raise_instead_of_silent_spectrum(csv, kwargs, match):
    """audit C-09: no more y_F = 0 / nan 'successful' results."""
    import io

    with pytest.raises(ValueError, match=match):
        run_microdosimetry(DeviceConfig(), mc_csv_path=io.BytesIO(csv), **kwargs)

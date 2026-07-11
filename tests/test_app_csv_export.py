"""Unit tests for the pure app.components.results.to_csv_bytes seam.

UI-06 coverage: CSV serialization for all three sim_types (cv, cce, field)
plus the commented metadata header and the unknown-sim_type error path. No
Streamlit runtime and no devsim build are required — to_csv_bytes is a pure
function operating on a hand-built SimResult fixture, imported directly
(mirrors tests/test_app_device_sidebar.py's pure-seam shape).
"""

from io import StringIO

import numpy as np
import pandas as pd
import pytest

from petringa import DeviceConfig, SimResult
from app.components.results import to_csv_bytes


def _decode(result: SimResult) -> str:
    return to_csv_bytes(result).decode("utf-8")


def test_cv_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.array([0.0, -5.0, -10.0])
    y = np.array([1e-9, 5e-10, 2e-10])
    one_over_c2 = np.array([1e17, 4e17, 25e17])
    depletion_widths = np.array([1e-4, 2e-4, 3e-4])
    result = SimResult(
        config=cfg,
        sim_type="cv",
        x=x,
        y=y,
        metadata={
            "one_over_C_squared": one_over_c2,
            "depletion_widths": depletion_widths,
        },
    )

    text = _decode(result)

    assert "# software_version: 5.0.0" in text
    device_line = [ln for ln in text.splitlines() if ln.startswith("# device:")]
    assert device_line, "missing '# device:' header line"
    assert "epi_thickness_um=" in device_line[0]

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert (
        non_comment_lines[0]
        == "bias_V,capacitance_F,one_over_C2_cm4_per_F2,depletion_width_cm"
    )

    # A data row matches the fixture values (first row: bias=0.0).
    first_data_row = non_comment_lines[1].split(",")
    assert float(first_data_row[0]) == pytest.approx(0.0)
    assert float(first_data_row[1]) == pytest.approx(1e-9)
    assert float(first_data_row[2]) == pytest.approx(1e17)
    assert float(first_data_row[3]) == pytest.approx(1e-4)

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == [
        "bias_V",
        "capacitance_F",
        "one_over_C2_cm4_per_F2",
        "depletion_width_cm",
    ]
    assert len(df) == 3


def test_cce_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.array([-10.0, -50.0, -100.0])
    y = np.array([0.5, 0.9, 1.0])
    I_collected = np.array([1e-6, 2e-6, 3e-6])
    I_generated = 3.2e-6  # scalar total generated current (confirmed via source read)
    result = SimResult(
        config=cfg,
        sim_type="cce",
        x=x,
        y=y,
        metadata={"I_collected": I_collected, "I_generated": I_generated},
    )

    text = _decode(result)

    assert "# software_version: 5.0.0" in text
    device_line = [ln for ln in text.splitlines() if ln.startswith("# device:")]
    assert device_line
    assert "epi_thickness_um=" in device_line[0]

    i_gen_lines = [
        ln for ln in text.splitlines() if ln.startswith("# I_generated_A_per_cm2:")
    ]
    assert i_gen_lines, "missing '# I_generated_A_per_cm2:' header line"
    assert "3.2e-06" in i_gen_lines[0] or "3.2e-6" in i_gen_lines[0]

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert non_comment_lines[0] == "bias_V,CCE,I_collected_A_per_cm2"
    assert "I_generated" not in non_comment_lines[0]

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == ["bias_V", "CCE", "I_collected_A_per_cm2"]
    assert "I_generated" not in df.columns
    assert len(df) == 3


def test_field_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.array([0.0, 5.0, 10.0])  # depth in um already
    y = np.array([1e5, 5e4, 0.0])
    potential = np.array([0.0, -100.0, -200.0])
    net_doping = np.array([1e19, 1e15, 1e15])
    result = SimResult(
        config=cfg,
        sim_type="field",
        x=x,
        y=y,
        metadata={"potential": potential, "net_doping": net_doping},
    )

    text = _decode(result)

    assert "# software_version: 5.0.0" in text
    device_line = [ln for ln in text.splitlines() if ln.startswith("# device:")]
    assert device_line
    assert "epi_thickness_um=" in device_line[0]

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert (
        non_comment_lines[0]
        == "depth_um,ElectricField_V_per_cm,Potential_V,NetDoping_cm-3"
    )

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == [
        "depth_um",
        "ElectricField_V_per_cm",
        "Potential_V",
        "NetDoping_cm-3",
    ]
    assert len(df) == 3


def test_unknown_sim_type_raises_value_error():
    cfg = DeviceConfig()
    result = SimResult(
        config=cfg,
        sim_type="damage",
        x=np.array([1.0]),
        y=np.array([0.5]),
        metadata={},
    )

    with pytest.raises(ValueError):
        to_csv_bytes(result)

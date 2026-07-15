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

from etna import DeviceConfig, SimResult
from app.components.results import (
    build_dark_current_figure,
    sweep_results_to_csv_bytes,
    to_csv_bytes,
)


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


def test_damage_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.array([1e13, 1e14, 1e15])
    y = np.array([0.99, 0.95, 0.80])
    result = SimResult(
        config=cfg,
        sim_type="damage",
        x=x,
        y=y,
        metadata={"V_bias": -20.0, "energy_MeV": 5.6},
    )

    text = _decode(result)

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert non_comment_lines[0] == "fluence_p_per_cm2,CCE"

    v_bias_lines = [ln for ln in text.splitlines() if ln.startswith("# V_bias:")]
    assert v_bias_lines, "missing '# V_bias:' header line"
    assert "-20.0" in v_bias_lines[0]

    warning_lines = [ln for ln in text.splitlines() if ln.startswith("# WARNING:")]
    assert warning_lines, "missing '# WARNING:' header line"
    assert "data-blocked" in warning_lines[0]

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == ["fluence_p_per_cm2", "CCE"]
    assert len(df) == 3


def test_dark_current_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.array([250.0, 325.0, 400.0])
    y = np.array([1e-13, 5e-12, 2e-10])
    result = SimResult(
        config=cfg,
        sim_type="dark_current",
        x=x,
        y=y,
        metadata={
            "I_SRH": np.array([1e-13, 5e-12, 2e-10]),
            "I_TAT": np.array([-1e-14, -3e-13, -5e-12]),
            "I_SRV": np.array([0.0, 0.0, 0.0]),
        },
    )

    text = _decode(result)

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert non_comment_lines[0] == "T_K,I_total_A,I_SRH_A,I_TAT_A,I_SRV_A"

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == ["T_K", "I_total_A", "I_SRH_A", "I_TAT_A", "I_SRV_A"]
    assert len(df) == 3


def test_microdosimetry_csv_columns_and_header():
    cfg = DeviceConfig()
    x = np.geomspace(0.01, 9000, 5)
    y = np.array([0.1, 0.5, 1.2, 0.8, 0.05])
    result = SimResult(
        config=cfg,
        sim_type="microdosimetry",
        x=x,
        y=y,
        metadata={"y_F": 17.23, "y_D": 53.22, "l_bar_um": 20.0},
    )

    text = _decode(result)

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert non_comment_lines[0] == "y_keV_per_um,y_times_d_y"

    y_f_lines = [ln for ln in text.splitlines() if ln.startswith("# y_F_keV_per_um:")]
    assert y_f_lines, "missing '# y_F_keV_per_um:' header line"
    assert "17.23" in y_f_lines[0]

    y_d_lines = [ln for ln in text.splitlines() if ln.startswith("# y_D_keV_per_um:")]
    assert y_d_lines, "missing '# y_D_keV_per_um:' header line"
    assert "53.22" in y_d_lines[0]

    l_bar_lines = [ln for ln in text.splitlines() if ln.startswith("# l_bar_um:")]
    assert l_bar_lines, "missing '# l_bar_um:' header line"
    assert "20.0" in l_bar_lines[0]

    df = pd.read_csv(StringIO(text), comment="#")
    assert list(df.columns) == ["y_keV_per_um", "y_times_d_y"]
    assert len(df) == 5


def test_sweep_results_to_csv_bytes_shape():
    cfg = DeviceConfig()
    results = [
        SimResult(
            config=cfg,
            sim_type="cce",
            x=np.arange(3, dtype=float),
            y=np.arange(3, dtype=float),
            metadata={},
        ),
        SimResult(
            config=cfg,
            sim_type="cce",
            x=np.arange(4, dtype=float),
            y=np.arange(4, dtype=float),
            metadata={},
        ),
        SimResult(
            config=cfg,
            sim_type="cce",
            x=np.arange(5, dtype=float),
            y=np.arange(5, dtype=float),
            metadata={},
        ),
    ]
    values = [10.0, 15.0, 20.0]

    text = sweep_results_to_csv_bytes(results, "epi_thickness_um", values).decode(
        "utf-8"
    )

    swept_lines = [ln for ln in text.splitlines() if ln.startswith("# swept_values:")]
    assert swept_lines, "missing '# swept_values:' header line"

    version_lines = [
        ln for ln in text.splitlines() if ln.startswith("# software_version:")
    ]
    assert version_lines, "missing '# software_version:' header line"

    non_comment_lines = [ln for ln in text.splitlines() if not ln.startswith("#")]
    assert non_comment_lines[0] == "epi_thickness_um,x,y"

    df = pd.read_csv(StringIO(text), comment="#")
    assert "epi_thickness_um" in df.columns
    # One CSV concatenating all three runs: 3 + 4 + 5 = 12 data rows.
    assert len(df) == sum(len(r.x) for r in results) == 12


def test_unknown_sim_type_raises_value_error():
    cfg = DeviceConfig()
    result = SimResult(
        config=cfg,
        sim_type="not_a_real_sim_type",
        x=np.array([1.0]),
        y=np.array([0.5]),
        metadata={},
    )

    with pytest.raises(ValueError):
        to_csv_bytes(result)


def test_build_dark_current_figure_guards_zero_and_negative():
    cfg = DeviceConfig()
    result = SimResult(
        config=cfg,
        sim_type="dark_current",
        x=np.array([250.0, 325.0, 400.0]),
        y=np.array([1e-13, 5e-12, 2e-10]),
        metadata={
            "I_SRH": np.array([1e-13, 5e-12, 2e-10]),
            "I_TAT": np.array([-1e-14, -3e-13, -5e-12]),
            "I_SRV": np.array([0.0, 0.0, 0.0]),
        },
    )

    fig = build_dark_current_figure(result)

    assert len(fig.data) == 3
    assert not any(t.name == "SRV (surface)" for t in fig.data)

    tat_trace = next(t for t in fig.data if t.name == "TAT (effective)")
    assert all(v >= 0 for v in tat_trace.y)

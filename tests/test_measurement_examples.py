"""SYNTHETIC example measurement files (examples/measurements/) and the plausibility /
unit checks of the measurement loaders (etna.measurements.plausibility)."""

import importlib.util
import pathlib
import warnings

import numpy as np
import pytest

from etna import measurements as M
from etna.measurements import MeasurementFileError
from etna.measurements.plausibility import (
    MeasurementPlausibilityWarning,
    plausibility_report,
)

EX = pathlib.Path(__file__).resolve().parents[1] / "examples" / "measurements"


def _module(name):
    spec = importlib.util.spec_from_file_location(f"_ex_{name}", EX / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def test_committed_examples_are_reproducible():
    built = _module("make_synthetic").build()
    for name, text in built.items():
        assert (EX / "data" / name).read_text() == text, name
        assert "SYNTHETIC" in text and "NOT a measurement" in text


def test_examples_round_trip():
    with warnings.catch_warnings():
        warnings.simplefilter("error", MeasurementPlausibilityWarning)
        res = _module("roundtrip").run()
    assert res["cv"]["V_min"] == -200.0 and res["cv"]["n"] == 101
    assert res["cv"]["rms_W_um"] < 1e-6
    assert res["sims"]["rms_log10_vs_truth"] < 1e-3
    assert res["sims"]["junction_source"] == "file metadata"
    assert res["sims"]["cv_rms_W_um_rebuilt"] < 1e-3
    assert len(res["it"]) == 2
    for d in res["it"].values():
        assert abs(d["delta_eV"]) < 1e-6


def test_examples_have_no_plausibility_warnings():
    data = EX / "data"
    for obj in (
        M.load_cv_csv(data / "cv_synthetic.csv"),
        M.load_sims_csv(data / "sims_synthetic.csv"),
        M.load_it_csv(data / "it_synthetic.csv"),
    ):
        assert plausibility_report(obj) == []


def test_realistic_sims_floor_censors_the_epi():
    text = (EX / "data" / "sims_synthetic.csv").read_text()
    text = text.replace("detection_limit_N_cm3: 1e+13", "detection_limit_N_cm3: 5e+14")
    lines = [
        (
            ln
            if ln.startswith("#") or ln.startswith("depth")
            else ",".join(
                [
                    ln.split(",")[0],
                    "0" if 0 < float(ln.split(",")[1]) < 5e14 else ln.split(",")[1],
                    ln.split(",")[2],
                ]
            )
        )
        for ln in text.splitlines()
    ]
    prof = M.load_sims_csv("\n".join(lines) + "\n")
    from etna.api.device import measured_cv_tabulated_config

    table = M.resample_profile(prof, measured_cv_tabulated_config(), method="interp")
    assert table.censored.sum() > 10
    with pytest.raises(ValueError, match="detection limit"):
        table.to_device_config(measured_cv_tabulated_config())


# ---------------------------------------------------------------- unit mix-ups
CV_HEAD = "# area_cm2: 0.04\nV_V,C_pF\n"


def _cv(rows, head=CV_HEAD):
    return head + "\n".join(f"{v},{c}" for v, c in rows) + "\n"


def test_farads_in_a_pf_column_is_refused_with_a_hint():
    rows = [(0, 2.0e-10), (-10, 6e-11), (-50, 3.4e-11)]
    with pytest.raises(MeasurementFileError, match="F in a pF column"):
        M.load_cv_csv(_cv(rows))


def test_pf_in_a_farad_column_is_refused_with_a_hint():
    rows = [(0, 207.0), (-10, 60.0), (-50, 34.0)]
    with pytest.raises(MeasurementFileError, match="pF in an F column"):
        M.load_cv_csv(_cv(rows, "# area_cm2: 0.04\nV_V,C_F\n"))


def test_impossible_area_is_refused():
    rows = [(0, 207.0), (-10, 60.0), (-50, 34.0)]
    with pytest.raises(MeasurementFileError, match="area"):
        M.load_cv_csv(_cv(rows, "# area_cm2: 400\nV_V,C_pF\n"))


def test_cv_capacitance_rising_with_reverse_bias_warns():
    rows = [(0, 207.0), (-10, 60.0), (-20, 80.0), (-50, 34.0)]
    with pytest.warns(MeasurementPlausibilityWarning, match="rises"):
        M.load_cv_csv(_cv(rows))


def test_sims_per_m3_is_refused():
    text = "depth_um,N_cm3\n0.5,3e19\n1.0,5e20\n2.0,1e23\n"
    with pytest.raises(MeasurementFileError, match="atomic density"):
        M.load_sims_csv(text)


def test_sims_log10_values_are_refused():
    text = "depth_um,N_cm3\n0.5,14.2\n1.0,14.5\n2.0,18.0\n"
    with pytest.raises(MeasurementFileError, match="log10"):
        M.load_sims_csv(text)


def test_sims_nm_in_an_um_column_is_refused():
    text = "depth_um,N_cm3\n100,1e14\n5000,1e14\n12000,1e18\n"
    with pytest.raises(MeasurementFileError, match="depth_nm"):
        M.load_sims_csv(text)


def test_it_degc_in_a_kelvin_column_warns():
    text = "V_V,T_K,I_pA\n-10,25,1\n-10,50,10\n-10,75,100\n"
    with pytest.warns(MeasurementPlausibilityWarning, match="degC"):
        M.load_it_csv(text)


def test_it_picoamps_in_an_ampere_column_is_refused():
    text = "V_V,T_K,I_A\n-10,300,4\n-10,320,240\n-10,340,9000\n"
    with pytest.raises(MeasurementFileError, match="not a dark current"):
        M.load_it_csv(text)


def test_plausibility_report_rejects_other_objects():
    with pytest.raises(TypeError):
        plausibility_report(np.zeros(3))

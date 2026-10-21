"""App alignment with the current library (2026-10-21 audit).

Opt-in incomplete ionization in the sidebar, the n+-aware cathode text, mixed-solver
badges, the Radiation Damage energy / kappa / trap texts, the Dark Current solver label
taken from the results, profile-aware Batch Sweep fields, labelled axes for every
generic plot, the microdosimetry chord model, the 'About this page' box on every page and
no stale physics text in the app sources. Facades are monkeypatched where a solve is not
the point; the rest runs for real on the batched engine (< 1 s).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np
import pytest
from streamlit.testing.v1 import AppTest

import etna
from etna import DeviceConfig, SimResult

APP = Path(__file__).resolve().parents[1] / "app"
TIMEOUT = 240


def _page_script(module: str, config_expr: str = "DeviceConfig()") -> str:
    return (
        "import streamlit as st, etna\n"
        "from etna import DeviceConfig\n"
        "if 'device_config' not in st.session_state:\n"
        f"    st.session_state['device_config'] = {config_expr}\n"
        f"from app.workflows.{module} import render\n"
        "render()\n"
    )


def _app(module: str, config_expr: str = "DeviceConfig()") -> AppTest:
    return AppTest.from_string(
        _page_script(module, config_expr), default_timeout=TIMEOUT
    )


# ------------------------------------------------------------------ sidebar


def test_assemble_config_passes_incomplete_ionization():
    from app.components.device_sidebar import assemble_config

    base = {
        "epi_thickness_um": 10.0,
        "substrate_thickness_um": DeviceConfig.substrate_thickness_um,
        "half_width_um": None,
        "N_A": 1e19,
        "doping_profile": "graded",
        "N_D": None,
        "N_D_junction": DeviceConfig.N_D_junction,
        "N_D_bulk": DeviceConfig.N_D_bulk,
        "L_transition_um": DeviceConfig.L_transition_um,
        "T": 300.0,
        "area_cm2": DeviceConfig.area_cm2,
    }
    assert assemble_config(base).incomplete_ionization is False
    assert assemble_config(
        {**base, "incomplete_ionization": True}
    ).incomplete_ionization


def test_sidebar_incomplete_ionization_toggle():
    at = AppTest.from_file(str(APP / "main.py"), default_timeout=TIMEOUT)
    at.run()
    box = at.sidebar.checkbox(key="cfg_incomplete_ionization")
    assert box.value is False and box.help
    assert at.session_state["device_config"].incomplete_ionization is False
    box.check()
    at.run()
    assert at.exception == []
    assert at.session_state["device_config"].incomplete_ionization is True
    # the sidebar says what it changes
    assert any(
        "Incomplete ionization on" in (c.value or "") for c in at.sidebar.caption
    )


def test_sidebar_labels_name_the_p_plus_anode():
    at = AppTest.from_file(str(APP / "main.py"), default_timeout=TIMEOUT)
    at.run()
    labels = [n.label for n in at.sidebar.number_input]
    assert "p+ anode thickness (µm)" in labels and "N_A p+ anode (cm⁻³)" in labels
    assert not any("substrate" in lab.lower() for lab in labels)
    # every sidebar number input explains itself
    assert all(n.help for n in at.sidebar.number_input), [
        n.label for n in at.sidebar.number_input if not n.help
    ]


# ------------------------------------------------------------------ page_utils


def _rt_result(cfg):
    return SimResult(
        config=cfg,
        sim_type="cv",
        x=np.array([0.0, -10.0]),
        y=np.array([1e-10, 5e-11]),
        metadata={"reach_through_V": 3.1},
    )


def test_cathode_text_follows_the_n_plus_layer():
    from app.components.page_utils import reach_through_text, truncation_text

    text = reach_through_text(_rt_result(DeviceConfig()))
    assert "n⁺ layer under the epi" in text and "directly on the n⁻ epi" not in text
    legacy = _rt_result(etna.legacy_config())
    assert "directly on the n⁻ epi" in reach_through_text(legacy)
    legacy.metadata.update(truncated=True, requested_v_stop=-100.0)
    assert "directly on the n⁻ epi" in truncation_text(legacy)


def test_solver_badge_names_mixed_solvers():
    from app.components.page_utils import solver_badge_text

    text = solver_badge_text({"solver": "batched + devsim"})
    assert "batched 1D engine (numba) + devsim (reference finite-volume solver)" in text
    assert "points from different solvers" in text


def test_trap_model_badge_is_honest():
    from app.components.page_utils import solver_badge_text

    text = solver_badge_text({"solver": "batched+traps"})
    assert "charged deep traps" in text and "no devsim reference" in text
    assert "not cross-checked against the batched engine" not in text


# ------------------------------------------------------------------ radiation damage


def test_radiation_energy_kappa_and_trap_texts(monkeypatch):
    from etna.api.simulation import compensation_fluence, trap_compensation_fluence
    from etna.core.radiation_damage import get_hardness_factor

    monkeypatch.setattr(etna, "run_radiation_damage", lambda cfg, **k: None)
    at = _app("radiation_damage")
    at.run()
    energy = at.number_input(key="rad_proton_energy")
    energy.set_value(5.6)
    at.run()
    assert at.exception == []
    k = get_hardness_factor(5.6)
    assert any(f"{k:.3g}" in (c.value or "") for c in at.caption)
    phi_c = compensation_fluence(DeviceConfig(), 5.6)
    assert any(f"{phi_c:.2e}" in (i.value or "") for i in at.info)
    at.checkbox(key="rad_traps").check()
    at.run()
    phi_t = trap_compensation_fluence(DeviceConfig(radiation_traps=True), 5.6)
    assert any(f"{phi_t:.2e}" in (i.value or "") for i in at.info)


def test_radiation_traps_disabled_with_incomplete_ionization(monkeypatch):
    monkeypatch.setattr(etna, "run_radiation_damage", lambda cfg, **k: None)
    at = _app("radiation_damage", "DeviceConfig(incomplete_ionization=True)")
    at.run()
    assert at.exception == []
    assert at.checkbox(key="rad_traps").disabled
    assert any("exclusive" in (c.value or "") for c in at.caption)


def test_refused_points_suggest_the_trap_model(monkeypatch):
    def _fake(cfg, fluences, **kwargs):
        reason = ["", "full donor compensation: outside model validity"]
        return SimResult(
            config=cfg,
            sim_type="damage",
            x=np.array([1e11, 1e13]),
            y=np.array([0.99, np.nan]),
            metadata={
                "V_bias": -20.0,
                "energy_MeV": 62.0,
                "valid": np.array([True, False]),
                "reason": reason,
            },
        )

    monkeypatch.setattr(etna, "run_radiation_damage", _fake)
    at = _app("radiation_damage")
    at.run()
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert any("Charged deep traps" in (w.value or "") for w in at.warning)


# ------------------------------------------------------------------ dark current


def test_dark_current_solver_label_from_results(monkeypatch):
    def _fake(cfg, v_start, v_stop, n_points, N_t=None, **kw):
        solver = "batched" if cfg.T < 450 else "devsim"
        return SimResult(
            config=cfg,
            sim_type="dark_current",
            x=np.array([v_stop]),
            y=np.array([1e-11]),
            metadata={"solver": solver},
        )

    monkeypatch.setattr(etna, "run_dark_current", _fake)
    at = _app("dark_current")
    at.run()
    at.number_input(key="dc_t_max").set_value(500.0)
    at.number_input(key="dc_n_temps").set_value(2)
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert at.session_state["dark_current_result"].metadata["solver"] == (
        "batched + devsim"
    )


# ------------------------------------------------------------------ batch sweep


def test_sweepable_fields_follow_the_profile():
    from app.workflows.batch_sweep import DEFAULT_VALUES, sweepable_fields

    graded = sweepable_fields(DeviceConfig())
    assert {"N_D_junction", "N_D_bulk", "L_transition_um"} <= set(graded)
    assert {"buffer_thickness_um", "buffer_end_doping"} <= set(graded)
    assert "N_D" not in graded
    uniform = sweepable_fields(DeviceConfig(doping_profile="uniform", N_D=1e15))
    assert "N_D" in uniform and "N_D_junction" not in uniform
    table = sweepable_fields(etna.measured_cv_tabulated_config())
    assert not {"N_D_junction", "N_D_bulk", "L_transition_um", "N_D"} & set(table)
    for cfg in (DeviceConfig(), DeviceConfig(doping_profile="uniform", N_D=1e15)):
        for f in sweepable_fields(cfg):
            values = [float(v) for v in DEFAULT_VALUES[f].split(",")]
            for v in values:  # every default value is a valid device
                import dataclasses

                dataclasses.replace(cfg, **{f: v})


def test_batch_sweep_uniform_device_offers_n_d():
    at = _app("batch_sweep", "DeviceConfig(doping_profile='uniform', N_D=1e15)")
    at.run()
    assert at.exception == []
    options = at.selectbox(key="sweep_param").options
    assert "N_D uniform (cm⁻³)" in options
    assert "N_D junction (cm⁻³)" not in options


# ------------------------------------------------------------------ plots


@pytest.mark.parametrize("sim_type", ["temperature", "flash", "transient", "sweep"])
def test_generic_results_have_labelled_axes(sim_type):
    from app.components.results import build_result_figure

    res = SimResult(
        config=DeviceConfig(),
        sim_type=sim_type,
        x=np.array([1.0, 2.0]),
        y=np.array([0.9, 0.95]),
        metadata={},
    )
    fig = build_result_figure(res, "t")
    xt, yt = fig.layout.xaxis.title.text, fig.layout.yaxis.title.text
    assert xt not in ("x", None) and yt not in ("y", None)
    assert "(" in xt or "CCE" in yt  # a unit or a named dimensionless quantity


def test_run_all_curves_use_the_dedicated_builders():
    from app.workflows.run_all import _curve

    report = etna.run_all(DeviceConfig(), ("cv", "temperature"))
    titles = {
        s.name: _curve(s).layout.xaxis.title.text
        for s in report.sections
        if s.status == "ok"
    }
    assert titles == {"cv": "Voltage (V)", "temperature": "Temperature (K)"}


def test_figures_depletion_view_plots_widths():
    from app.workflows.figures import _plotly

    res = etna.run_cv(DeviceConfig(), 0.0, -10.0, 6)
    fig = _plotly(res, "depletion")
    assert fig.layout.yaxis.title.text == "Depletion Width (µm)"
    assert len(fig.data) == 2
    assert _plotly(res, "map_electrons") is None


# ------------------------------------------------------------------ microdosimetry


def test_microdosimetry_passes_the_chord_model(monkeypatch):
    seen = {}

    def _fake(cfg, mc_csv_path, sv_thickness_um, sv_width_um, chord_model="cauchy"):
        seen["chord_model"] = chord_model
        return SimResult(
            config=cfg,
            sim_type="microdosimetry",
            x=np.array([1.0, 10.0]),
            y=np.array([0.1, 0.2]),
            metadata={"y_F": 1.0, "y_D": 2.0, "l_bar_um": 20.0},
        )

    monkeypatch.setattr(etna, "run_microdosimetry", _fake)
    at = _app("microdosimetry")
    at.run()
    data = Path(__file__).resolve().parents[1] / "data" / "synthetic_mc_events.csv"
    at.file_uploader[0].upload("events.csv", data.read_bytes(), "text/csv")
    at.selectbox(key="micro_chord").set_value("slab")
    at.run()
    at.button[0].click()
    at.run()
    assert at.exception == []
    assert seen["chord_model"] == "slab"
    assert any("l̄ = 20.00 µm" in (c.value or "") for c in at.caption)


# ------------------------------------------------------------------ every page


PAGES = [
    "cv",
    "cce",
    "field_map",
    "radiation_damage",
    "dark_current",
    "microdosimetry",
    "batch_sweep",
    "run_all",
    "export",
    "figures",
    "measurements",
]


@pytest.mark.parametrize("module", PAGES)
def test_every_page_explains_itself(module):
    at = _app(module)
    at.run()
    assert at.exception == []
    labels = [e.label for e in at.expander]
    assert any(lab.startswith("About this page") for lab in labels), labels


def test_no_stale_physics_text_in_the_app():
    stale = (
        "no n⁺ substrate",  # the default stack has the n+ layer
        "placeholder kappa",
        "kappa > 1",
        "κ > 1",
        "3-point n+ calibration",  # the previous preset is a 2-point C-V fit
        "set the half-width to none",
        'xaxis_title="x"',
    )
    for path in APP.rglob("*.py"):
        text = path.read_text()
        for s in stale:
            assert s not in text, (path.name, s)

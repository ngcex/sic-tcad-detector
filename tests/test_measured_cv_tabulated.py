"""The OPT-IN tabulated donor profile (WP-E): the 19-knot log-linear inversion of the
digitised measured C-V as a DeviceConfig.donor_table, the devsim builders (1D/2D), the
batched 1D engine, the 2D engine gate, the preset measured_cv_tabulated_config(), the
full-depletion definitions, the exports and the sidebar.

Not the default: DEFAULT_STACK / DeviceConfig() are unchanged (checked first). Numbers:
lab/tabulated_profile/notes/LOG.md (1D parity dW 1.4e-15, dpsi 1.4e-14 V, dCCE 2.8e-11 on
the preset; regression vs the digitised curve knee rms 0.015 um, plateau rms 0.0075 um).
"""

import dataclasses
import json
import uuid

import devsim
import numpy as np

from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM
import pytest

import etna
from etna.api import simulation as sim
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    measured_cv_config,
    measured_cv_tabulated_config,
)
from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.batched.material import EPS
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.device import (
    DEFAULT_STACK,
    MEASURED_CV_TABLE,
    MEASURED_CV_TABULATED_STACK,
    donor_table_spec,
    mobility_doping,
    table_min_doping,
    table_sensitive_end_cm,
    table_sensitive_end_doping,
    tabulated_donor_profile,
)
from etna.core.drift_diffusion import create_dd_device

S = MEASURED_CV_TABULATED_STACK
SUB, EPI, NP = S["substrate_thickness_cm"], 10e-4, S["nplus_thickness_cm"]
SPEC = donor_table_spec(EPI, MEASURED_CV_TABLE)


def _kw(**over):
    kw = device_kwargs(measured_cv_tabulated_config())
    kw.update(over)
    return kw


# ------------------------------------------------------------- opt-in: the default is intact
def test_default_stack_unchanged_and_preset_record():
    assert "donor_table" not in DEFAULT_STACK
    assert DeviceConfig().donor_table is None and DeviceConfig() == measured_cv_config()
    assert "donor_table" not in device_kwargs(DeviceConfig())
    c = measured_cv_tabulated_config()
    assert c != DeviceConfig()
    assert c.donor_table == MEASURED_CV_TABLE and len(c.donor_table) == 19
    assert c.doping_profile == "graded"
    assert c.buffer_thickness_um is None and c.buffer_end_doping is None
    assert (c.substrate_thickness_um, c.nplus_thickness_um, c.N_D_nplus) == (
        0.3,
        1.0,
        1e19,
    )
    assert c.N_t == S["N_t"] == 1.356120e12 and c.E_t == 0.67
    # informative doping fields = the table's junction value and minimum
    assert (
        c.N_D_junction == MEASURED_CV_TABLE[0][1] == pytest.approx(5.927e15, rel=1e-3)
    )
    assert c.N_D_bulk == table_min_doping(SPEC) == pytest.approx(1.295e13, rel=1e-3)
    # the knots are the lab inversion (knots.json, lam 1, no fringe), bit for bit
    with open("lab/cv_gap/notes/knots.json") as f:
        lab = json.load(f)["knots lam 1.0 Qf None"]
    assert [n for _, n in MEASURED_CV_TABLE] == [10.0**l for l in lab["logN_knots"]]
    assert [d for d, _ in MEASURED_CV_TABLE] == [
        0.0, 0.25, 0.5, 0.75, 1.0, 1.5, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0, 8.5, 9.0,
        9.25, 9.5, 9.75, 10.0,
    ]  # fmt: skip
    assert etna.MEASURED_CV_TABULATED_PRESET["donor_table"] is MEASURED_CV_TABLE
    assert measured_cv_tabulated_config(T=310.0).T == 310.0
    # device_kwargs carries the table and explicit Nones for the transition
    kw = device_kwargs(c)
    assert kw["donor_table"] == MEASURED_CV_TABLE
    assert kw["buffer_thickness_cm"] is None and kw["buffer_end_doping"] is None


def test_config_validation():
    tab = ((0.0, 1e15), (5.0, 5e13), (10.0, 1e16))
    # lists are normalised to a tuple of float pairs; the default transition is dropped
    c = DeviceConfig(donor_table=[[0, 1e15], [5, 5e13], [10, 1e16]])
    assert c.donor_table == tab and c.buffer_thickness_um is None
    with pytest.raises(ValueError, match="graded"):
        DeviceConfig(doping_profile="uniform", N_D=1e15, donor_table=tab)
    with pytest.raises(ValueError, match="exclusive"):
        DeviceConfig(donor_table=tab, buffer_thickness_um=1.0, buffer_end_doping=1e16)
    with pytest.raises(ValueError, match="cover the epi"):
        DeviceConfig(donor_table=tab, epi_thickness_um=12.0)
    with pytest.raises(ValueError, match="increasing"):
        DeviceConfig(donor_table=((0.0, 1e15), (5.0, 5e13), (4.0, 1e16), (10.0, 1e16)))
    with pytest.raises(ValueError, match="start at depth 0"):
        DeviceConfig(donor_table=((1.0, 1e15), (10.0, 1e16)))
    with pytest.raises(ValueError, match="finite"):
        DeviceConfig(donor_table=((0.0, 1e15), (10.0, -1.0)))
    with pytest.raises(ValueError, match="two"):
        DeviceConfig(donor_table=((0.0, 1e15),))
    with pytest.raises(ValueError, match="below N_A"):
        DeviceConfig(donor_table=((0.0, 2e19), (10.0, 1e16)))
    # the core builders refuse the same combinations
    with pytest.raises(ValueError, match="exclusive"):
        create_dd_device(
            device_name="never",
            donor_table=tab,
            buffer_thickness_cm=1e-4,
            buffer_end_doping=1e16,
        )
    with pytest.raises(ValueError, match="graded"):
        create_dd_device(
            device_name="never", doping_profile="uniform", N_D=1e15, donor_table=tab
        )


def test_table_helpers():
    x = np.array([0.0, 0.3e-4, 0.3e-4 + 2.5e-4, 0.3e-4 + 10e-4, 0.3e-4 + 11e-4])
    d = tabulated_donor_profile(x, 0.3e-4, SPEC)
    assert d[0] == 0.0 and d[1] == MEASURED_CV_TABLE[0][1]
    # log-linear between the 2 um (14.234) and 3 um (13.701) knots
    assert d[2] == pytest.approx(10 ** ((14.23424292213302 + 13.701401160537118) / 2))
    assert d[3] == d[4] == MEASURED_CV_TABLE[-1][1]  # last value held beyond the table
    assert tabulated_donor_profile(x, 0.3e-4, SPEC, n_side=False)[0] == d[1]
    # sensitive-epi end: one decade above the minimum, on the final rise (9.21 um)
    assert table_sensitive_end_cm(SPEC) * 1e4 == pytest.approx(9.207, abs=0.002)
    assert table_sensitive_end_doping(SPEC) == pytest.approx(
        10 * table_min_doping(SPEC)
    )
    flat = donor_table_spec(10e-4, ((0.0, 1e14), (10.0, 1e14)))
    assert table_sensitive_end_cm(flat) == 10e-4  # never rises: the whole table
    assert mobility_doping(None, "graded", 3e13, SPEC) == table_min_doping(SPEC)


# ------------------------------------------------------- devsim doping == numpy twin (bits)
def test_devsim_donors_are_the_numpy_twin_bit_for_bit():
    name = f"tab_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **_kw())
    try:
        g = lambda n: np.array(
            devsim.get_node_model_values(device=name, region="sic", name=n)
        )
        x, d0 = g("x"), g("Donors")
        cv_sweep(info, np.array([0.0, -5.0, -30.0]))
        d1, net = g("Donors"), g("NetDoping")
    finally:
        devsim.delete_device(device=name)
    twin = tabulated_donor_profile(x, SUB, SPEC) + 1e19 * (x >= SUB + EPI)
    assert np.array_equal(d0, twin) and np.array_equal(d1, twin)
    assert np.array_equal(net[x > SUB], twin[x > SUB])
    assert info["donor_table"] == MEASURED_CV_TABLE
    assert info["buffer_thickness_cm"] is None and info["buffer_end_doping"] is None
    assert info["N_D_mobility"] == table_min_doping(SPEC)


# ------------------------------------------------------------ parity engine vs devsim (1D)
BIASES = np.array([0.0, -2.0, -10.0, -40.0, -100.0])


def test_engine_matches_devsim_with_table():
    kw = _kw()
    name = f"tabp_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **kw)
    try:
        ref = cv_sweep(info, BIASES)
        psi_ref = np.array(
            devsim.get_node_model_values(device=name, region="sic", name="Potential")
        )
    finally:
        devsim.delete_device(device=name)
    assert not ref["truncated"]
    x = node_positions(EPI, SUB, NP)
    pb = make_batch(
        x,
        SUB,
        1e19,
        1e15,
        1e14,
        1e-4,
        nplus=(SUB + EPI, 1e19),
        buffer=None,
        donor_table=SPEC,
    )
    res = reverse_sweep(pb, BIASES, SUB, EPI, keep_psi=True)
    assert res.converged.all()
    np.testing.assert_allclose(res.W[:, 0], ref["depletion_widths"], rtol=1e-10)
    np.testing.assert_allclose(res.psi[-1, 0], psi_ref, atol=1e-10)
    c = cce_vs_bias(BIASES, device_kwargs=kw)
    G = alpha_generation(x, SUB, ALPHA_RANGE_SIC_CM, 1e18, epi_end=SUB + EPI)
    ce = cce_sweep(pb, BIASES, SUB, G)
    assert ce.converged.all()
    np.testing.assert_allclose(ce.cce[:, 0], c["cce_values"], atol=1e-10)
    with pytest.raises(ValueError, match="exclusive"):
        make_batch(
            x, SUB, 1e19, 1e15, 1e14, 1e-4,
            nplus=(SUB + EPI, 1e19), buffer=(SUB + EPI, 1e-4, 1e16), donor_table=SPEC,
        )  # fmt: skip


# ------------------------------------------------------------------- engine boxes
def test_preset_inside_every_engine_domain():
    from etna.core.engine2d import (
        ENGINE2D_TABLE_BOX,
        engine2d_dark_unsupported_reason,
        engine2d_unsupported_reason,
    )

    c = measured_cv_tabulated_config()
    assert sim.validated_domain_reason(c) is None
    kw = device_kwargs(c)
    kw.pop("N_D"), kw.pop("epi_thickness_cm")
    nd = kw.pop("N_D_bulk")
    for hw in (25.0, 50.0):
        for vb in (20.0, 40.0, 60.0):
            assert engine2d_unsupported_reason(hw, vb, 10e-4, nd, 10, **kw) is None
    assert engine2d_dark_unsupported_reason(50.0, 50.0, N_D_bulk=nd, **kw) is None
    # the table's N_D range is the box (parametric N_D_bulk / N_D_junction / L not checked:
    # the preset's informative N_D_bulk 1.3e13 is below the parametric box's 2e13)
    assert nd < sim._ENGINE_BOX["N_D_bulk"][0]
    assert sim._ENGINE_BOX["donor_table_N_D"] == ENGINE2D_TABLE_BOX["donor_table_N_D"]
    low = tuple((d, n * 10.0**-0.5) for d, n in MEASURED_CV_TABLE)  # min 4.1e12
    c_low = measured_cv_tabulated_config(donor_table=low)
    assert "donor_table" in sim.validated_domain_reason(c_low)
    kl = device_kwargs(c_low)
    kl.pop("N_D"), kl.pop("epi_thickness_cm"), kl.pop("N_D_bulk")
    assert "donor_table" in engine2d_unsupported_reason(25.0, 20.0, 10e-4, nd, 3, **kl)
    assert "donor_table" in engine2d_dark_unsupported_reason(50.0, 50.0, **kl)
    # a table with an explicit transition, or a uniform profile, is refused by the gate
    assert "exclusive" in engine2d_unsupported_reason(
        25.0, 20.0, 10e-4, nd, 3,
        **{**kw, "buffer_thickness_cm": 1e-4, "buffer_end_doping": 1e16},
    )  # fmt: skip
    assert "graded" in engine2d_unsupported_reason(
        25.0, 20.0, 10e-4, nd, 3, **{**kw, "doping_profile": "uniform"}
    )
    # the gate resolves the builder's default transition as "none" with a table
    kw_nobuf = {k: v for k, v in kw.items() if not k.startswith("buffer")}
    assert engine2d_unsupported_reason(25.0, 20.0, 10e-4, nd, 3, **kw_nobuf) is None


# ------------------------------------------- the digitised measured C-V (fit reproduction)
def _digitised():
    root = "lab/cv_gap/notes/"
    try:
        band = np.genfromtxt(root + "photons_w.csv", delimiter=",", names=True)
        mk = np.genfromtxt(root + "photons_w_markers.csv", delimiter=",", names=True)
    except OSError:  # pragma: no cover
        pytest.skip("digitised C-V not available")
    plateau = [
        (float(v), float(w))
        for v, w, wt in zip(band["V"], band["W_um"], band["weight"])
        if np.isfinite(w)
        and wt >= 10
        and v <= -5.0
        and v not in (-10.0, -30.0, -9.75, -10.25, -10.5, -9.5, -29.75, -30.25)
    ]
    knee = [
        (round(float(v) * 4) / 4, float(w))
        for v, w in zip(mk["V"], mk["W_um"])
        if -3.3 < v < -0.1
    ]
    return plateau, knee


def test_preset_reproduces_the_digitised_cv():
    """GATES: knee (13 markers, -3.25..-0.25 V) rms <= 0.05 um and plateau (-5..-59 V) rms
    <= 0.02 um through run_cv (eps*A/C of the small-signal C, A = 4 mm2); every knee marker
    within its digitisation sigma (0.03 um (+) 0.04 V x |dW/dV|). The parametric default
    leaves 0.24 um rms in the knee (tests/test_measured_cv_preset.py)."""
    cfg = measured_cv_tabulated_config(area_cm2=0.04)
    V = np.round(np.arange(0.0, -60.001, -0.25), 2)
    r = sim.run_cv(cfg, 0.0, -60.0, len(V))
    assert r.metadata["solver"] == "batched" and r.metadata["validated_domain"]
    W = EPS * cfg.area_cm2 / r.y * 1e4
    at = lambda v: W[np.argmin(np.abs(r.x - v))]
    # the lab inversion's own W (knots.json W_key). The lab model held the table's last
    # value (1.33e16) through the n+ region; the preset adds the 1e19 n+ step at the epi
    # end as every stack does, which moves W(-58 V) by 2.4e-3 um (nothing elsewhere).
    np.testing.assert_allclose(
        [at(0), at(-2), at(-10), at(-30), at(-58)],
        [1.678536, 7.187552, 9.527993, 9.739940, 9.886451],
        atol=3e-3,
    )
    # labels 1.68 / 9.52 / 9.73
    assert abs(at(0) - 1.68) < 0.01 and abs(at(-10) - 9.52) < 0.01
    assert abs(at(-30) - 9.73) < 0.015
    plateau, knee = _digitised()
    res_p = np.array([w - at(v) for v, w in plateau])
    assert len(plateau) > 200 and np.sqrt(np.mean(res_p**2)) <= 0.02
    res_k = np.array([w - at(v) for v, w in knee])
    assert len(knee) == 13 and np.sqrt(np.mean(res_k**2)) <= 0.05
    slope = np.gradient(W, r.x)
    sig = np.sqrt(
        0.03**2
        + (0.04 * np.array([slope[np.argmin(np.abs(r.x - v))] for v, _ in knee])) ** 2
    )
    assert np.all(np.abs(res_k) <= sig)  # 13/13 within 1 sigma (parametric: 2/13)
    assert abs(at(-2) - 7.18) < 0.02 and abs(at(-5) - 9.27) < 0.01
    # reach-through by the decade rule (end of the low-doped epi 9.21 um, L_D 0.33 um)
    assert r.metadata["reach_through_V"] == pytest.approx(4.85, abs=0.051)


# ------------------------------------------------- full-depletion definitions for a table
def test_full_depletion_definitions_for_table():
    from etna.core.analytical import built_in_potential
    from etna.core.optimization import analytic_vfd_from_info
    from etna.export.labels import (
        VFD_DEFINITION,
        analytic_vfd,
        reach_through_definition,
        table_sensitive_end_um,
    )

    cfg = measured_cv_tabulated_config()
    assert table_sensitive_end_um(cfg) == pytest.approx(9.207, abs=0.002)
    v = analytic_vfd(cfg)
    # moment integral to 9.21 um: between the parametric 2.74 V (to 8.86 um) and the
    # engine's 4.85 V reach-through, consistent with the measured knee at -4..-5 V
    assert 2.74 < v < 4.85 and v == pytest.approx(3.373, abs=0.02)
    assert "decade" in reach_through_definition(cfg) and "TABULATED" in VFD_DEFINITION
    # the 2D sweep's gate (device_info based) and the export agree
    name = f"tabv_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **_kw())
    try:
        v_bi = built_in_potential(info["N_A_ionized"], cfg.N_D_junction, info["n_i"])
        assert analytic_vfd_from_info(info, v_bi) == pytest.approx(v, rel=1e-6)
    finally:
        devsim.delete_device(device=name)
    # the compensation fluence uses the table minimum
    assert sim.compensation_fluence(cfg, 5.6) == pytest.approx(
        sim.compensation_fluence(DeviceConfig(), 5.6)
        * table_min_doping(SPEC)
        / DEFAULT_STACK["N_D_bulk"]
    )


# ------------------------------------------------------------------ exports and the app
def test_exports_label_the_tabulated_preset():
    from etna.export.cross_section import cross_section_svg, stack_layers
    from etna.export.labels import (
        TABULATED_CALIBRATION_NOTE,
        epi_donor_profile,
        is_measured_cv_tabulated_device,
        is_reference_device,
    )
    from etna.export.spec_sheet import build_spec, spec_to_json, spec_to_markdown

    cfg = measured_cv_tabulated_config(area_cm2=0.04)
    assert is_measured_cv_tabulated_device(cfg) and not is_reference_device(cfg)
    assert not is_measured_cv_tabulated_device(DeviceConfig())
    spec = build_spec(cfg, None)
    assert spec["measured_cv_tabulated_preset"] and not spec["reference_device"]
    assert spec["model"].startswith("opt-in tabulated measured-C-V profile")
    cal = spec["cv"]["calibration"]
    assert cal == TABULATED_CALIBRATION_NOTE
    for s in (
        "digitised Photons Fig. 4b",
        "chi2 32 / 211",
        "knee rms 0.014 um",
        "plateau rms 0.007 um",
        "inversion of ONE C-V curve",
        "NOT a SIMS",
        "+-0.2 dex",
        "first 0.5 um",
        "8.5-9.25 um",
    ):
        assert s in cal, s
    assert (
        "inversion" in spec["statements"][1] and "not a SIMS" in spec["statements"][1]
    )
    ed = spec["epi_doping"]
    assert (
        len(ed["donor_table"]) == 19
        and ed["donor_table"][0]["N_D_cm3"] == MEASURED_CV_TABLE[0][1]
    )
    assert ed["table_sensitive_epi_end_um"] == pytest.approx(9.207, abs=0.002)
    assert ed["buffer_thickness_um"] is None and "tabulated" in ed["formula"]
    assert "refitted on the tabulated profile" in spec["dark_current"]["note"]
    md = spec_to_markdown(spec)
    assert "| knot | depth below junction (um) | N_D (cm^-3) | log10 N_D |" in md
    assert "| 18 | 10 | 1.334e+16 | 16.125 |" in md
    json.loads(spec_to_json(spec))
    names = [L["name"] for L in stack_layers(cfg)]
    assert names == ["p+ anode", "n- epi", "n-/n+ transition", "n+ substrate"]
    layer = stack_layers(cfg)[2]
    assert layer["thickness_um"] == pytest.approx(10 - 9.207, abs=0.002)
    assert "not a SIMS" in layer["note"] and "tabulated" in layer["doping"]
    assert "transition" in cross_section_svg(cfg)
    assert float(epi_donor_profile(cfg, 5.0)) == table_min_doping(SPEC)
    assert float(epi_donor_profile(cfg, 10.0)) == pytest.approx(
        MEASURED_CV_TABLE[-1][1]
    )
    # a user table is labelled as such, without the preset's provenance
    other = DeviceConfig(donor_table=((0.0, 1e15), (5.0, 5e13), (10.0, 1e16)))
    sp = build_spec(other, None)
    assert sp["model"].startswith("tabulated epi profile (3 knots, user table)")
    assert not sp["measured_cv_tabulated_preset"]
    assert "provenance not known" in sp["epi_doping"]["table_note"]
    assert "reference device" in sp["cv"]["calibration"]


def test_sidebar_assembles_the_tabulated_preset():
    from app.components.device_sidebar import (
        FIELD_LABELS,
        assemble_config,
        humanize_config_error,
    )

    values = {
        "epi_thickness_um": 10.0,
        "substrate_thickness_um": 0.3,
        "half_width_um": None,
        "N_A": 1e19,
        "doping_profile": "graded",
        "N_D": None,
        "N_D_junction": 1e15,
        "N_D_bulk": 1e14,
        "L_transition_um": 1.0,
        "T": 300.0,
        "area_cm2": 0.04,
    }
    assert assemble_config(values).donor_table is None
    values.update(etna.MEASURED_CV_TABULATED_PRESET)
    assert assemble_config(values) == measured_cv_tabulated_config(area_cm2=0.04)
    assert "donor_table" in FIELD_LABELS
    assert "tabulated epi profile" in humanize_config_error("donor_table: bad")


def test_sidebar_toggle_applies_the_tabulated_preset():
    from app.components.device_sidebar import render_device_sidebar  # noqa: F401
    from streamlit.testing.v1 import AppTest

    at = AppTest.from_file("app/main.py", default_timeout=120)
    at.run()
    assert not at.sidebar.checkbox(key="cfg_tabulated_on").value
    assert at.session_state["device_config"].donor_table is None
    at.sidebar.checkbox(key="cfg_tabulated_on").check()
    at.run()
    assert at.session_state["device_config"] == dataclasses.replace(
        measured_cv_tabulated_config(), area_cm2=0.04
    )
    infos = " ".join(c.value for c in at.sidebar.info)
    assert "Tabulated measured-C-V profile active" in infos and "Not a SIMS" in infos
    texts = " ".join(c.value for c in at.sidebar.caption)
    assert "Inside the 1D engine/devsim cross-check domain" in texts
    # the previous calibrated model wins over the table when both are ticked (warned)
    at.sidebar.checkbox(key="cfg_nplus_cal_on").check()
    at.run()
    assert at.session_state["device_config"].donor_table is None
    assert any("Both presets" in w.value for w in at.sidebar.warning)


def test_run_all_sweep_on_the_engine_with_table():
    from etna.export.labels import analytic_vfd

    cfg = measured_cv_tabulated_config(area_cm2=0.04)
    rep = etna.run_all(cfg, sections=("sweep",))
    s = rep.get("sweep")
    assert s.status == "ok" and s.solver == "engine2d", (s.status, s.reason)
    assert "tabulated" in s.summary
    df = s.result.metadata["table"]
    assert df["V_fd"].iloc[0] == pytest.approx(analytic_vfd(cfg), rel=1e-3)
    assert int(df["is_valid"].sum()) == len(df) == 6


def test_batches_group_by_table_and_mix_with_the_default():
    tab = measured_cv_tabulated_config(area_cm2=0.04)
    other = measured_cv_tabulated_config(
        area_cm2=0.04, donor_table=tuple((d, n * 1.5) for d, n in MEASURED_CV_TABLE)
    )
    dflt = DeviceConfig(area_cm2=0.04)
    batch = sim.run_cv_batch([tab, dflt, other], 0.0, -40.0, 9)
    for r, c in zip(batch, (tab, dflt, other)):
        assert r.metadata["solver"] == "batched"
        single = sim.run_cv(c, 0.0, -40.0, 9)
        np.testing.assert_allclose(r.y, single.y, rtol=1e-12)
    assert not np.allclose(batch[0].y, batch[2].y, rtol=1e-3, atol=0.0)
    # 1D field map (batched engine and devsim) and the 2D dark state (engine2d) take the
    # table
    f = sim.run_field(tab, -30.0)
    assert f.metadata["solver"] == "batched" and np.max(f.y) > 5e4
    fd = sim.run_field(tab, -30.0, solver="devsim")
    np.testing.assert_allclose(f.y, fd.y, rtol=0, atol=5e-12 * fd.y.max())
    f2 = sim.run_field(dataclasses.replace(tab, half_width_um=50.0), -30.0)
    assert f2.metadata["solver"] == "engine2d"
    nd = np.asarray(f2.mesh.node_values["NetDoping"])
    assert nd.max() > 1e19 and np.isclose(
        np.sort(nd[(nd > 1e12) & (nd < 1e15)])[0], table_min_doping(SPEC), rtol=0.3
    )


# ---------------------------------------------------------------------------- slow
@pytest.mark.slow
def test_engine2d_matches_devsim_with_table():
    from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
    from etna.core.engine2d import lateral_cce_scan_2d

    kw = _kw()
    kw.pop("N_D"), kw.pop("doping_profile")
    e = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=3, **kw)
    assert e["device_info"]["donor_table"] == MEASURED_CV_TABLE
    name = f"tab2d_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=25.0, V_bias=20.0, device_name=name, **kw)
    try:
        y = np.array(devsim.get_node_model_values(device=name, region="sic", name="y"))
        d = np.array(
            devsim.get_node_model_values(device=name, region="sic", name="Donors")
        )
        scan = cce_lateral_scan(info, n_points=3)
    finally:
        devsim.delete_device(device=name)
    twin = tabulated_donor_profile(y, SUB, SPEC) + 1e19 * (y >= SUB + EPI)
    assert np.array_equal(d, twin)
    np.testing.assert_allclose(e["cce_values"], scan["cce_values"], rtol=0, atol=1e-10)


@pytest.mark.slow
def test_nt_reproduces_the_18pA_bound_on_the_tabulated_preset():
    from etna.core.dark_current import (
        AREA_DOSIMETER_CM2,
        create_dark_current_device,
        dark_current_sweep,
        setup_tat_model,
    )

    info = create_dark_current_device(
        N_t=S["N_t"], E_t=S["E_t"], device_name=f"tabnt_{uuid.uuid4().hex[:6]}", **_kw()
    )
    try:
        setup_tat_model(info, E_t=S["E_t"], N_t=S["N_t"])
        r = dark_current_sweep(info, np.array([-60.0]), area=AREA_DOSIMETER_CM2)
    finally:
        devsim.delete_device(device=info["device_name"])
    assert abs(r["I_total"][0]) == pytest.approx(18e-12, rel=1e-5)

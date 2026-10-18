"""The measured-C-V profile (graded n-/n+ transition, "buffer"): the DEFAULT since
2026-10-18 (it was an opt-in preset; measured_cv_config() is now an alias of DeviceConfig()).

Record: etna.core.device.DEFAULT_STACK = MEASURED_CV_STACK (lab/cv_gap fit on the full
digitised Photons Fig. 4b curve). Checks: the previous default (nplus_calibrated_config,
no buffer) is bit-identical to its closed formula; devsim and the
batched engine agree with the buffer (lab/buffer_profile/bench/buffer_parity1d.py: dW
1.9e-15, dpsi 1.4e-14 V, dCCE 1.2e-10 over 11 cases to -100 V; 2D: buffer_parity2d.py dpsi
7.1e-14 V, dCCE 1.4e-14); the preset reproduces the digitised plateau (rms 0.013 um) and
the three labels; exports and the app label it as a fit (fit input, not a validation).
"""

import uuid

import devsim
import numpy as np
import pytest

import etna
from etna.api import simulation as sim
from etna.api.device import (
    DeviceConfig,
    device_kwargs,
    measured_cv_config,
    nplus_calibrated_config,
)
from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.batched.material import EPS, graded_doping
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.device import (
    DEFAULT_STACK,
    MEASURED_CV_STACK,
    NPLUS_CALIBRATION,
    donor_profile,
)
from etna.core.drift_diffusion import create_dd_device

M = MEASURED_CV_STACK
SUB, EPI, NP = M["substrate_thickness_cm"], 10e-4, M["nplus_thickness_cm"]
BUFFER = (SUB + EPI, M["buffer_thickness_cm"], M["buffer_end_doping"])


def _preset_kwargs(**over):
    kw = device_kwargs(measured_cv_config())
    kw.update(over)
    return kw


# ------------------------------------------------- the default flip; previous default intact
def test_default_is_the_preset_and_previous_default_has_no_buffer():
    assert MEASURED_CV_STACK is DEFAULT_STACK
    assert DeviceConfig() == measured_cv_config()
    c = nplus_calibrated_config()
    assert c.buffer_thickness_um is None and c.buffer_end_doping is None
    kw = device_kwargs(c)
    # explicit Nones: the builders default to the stack's transition since the flip
    assert kw["buffer_thickness_cm"] is None and kw["buffer_end_doping"] is None
    P = NPLUS_CALIBRATION
    # the engine's doping without a buffer is the old closed formula, bit for bit
    x = node_positions(EPI, SUB, NP)
    d, _ = graded_doping(
        x,
        SUB,
        1e19,
        P["N_D_junction"],
        P["N_D_bulk"],
        P["L_transition"],
        nplus=(SUB + EPI, 1e19),
    )
    old = (
        P["N_D_bulk"]
        + (P["N_D_junction"] - P["N_D_bulk"])
        * np.exp(-np.maximum(x - SUB, 0.0) / P["L_transition"])
    ) * (x >= SUB) + 1e19 * (x - (SUB + EPI) >= 0)
    assert np.array_equal(d[0], old)


def test_preset_record_and_validation():
    c = measured_cv_config()
    assert (c.N_D_junction, c.N_D_bulk, c.L_transition_um) == (
        8.185880e15,
        3.382506e13,
        0.463863,
    )
    assert (c.buffer_thickness_um, c.buffer_end_doping) == (1.136634, 4.516421e16)
    assert (c.substrate_thickness_um, c.nplus_thickness_um, c.N_D_nplus) == (
        0.3,
        1.0,
        1e19,
    )
    assert c.N_t == M["N_t"] == 2.849024e12 and c.E_t == 0.67
    assert measured_cv_config(T=310.0).T == 310.0
    assert etna.MEASURED_CV_PRESET["buffer_end_doping"] == M["buffer_end_doping"]
    with pytest.raises(ValueError, match="together"):
        DeviceConfig(buffer_thickness_um=1.0, buffer_end_doping=None)
    with pytest.raises(ValueError, match="together"):
        DeviceConfig(buffer_thickness_um=None, buffer_end_doping=1e16)
    # the default transition (1.14 um) does not fit a 1 um epi: the message says how out
    with pytest.raises(ValueError, match="abrupt epi end"):
        DeviceConfig(epi_thickness_um=1.0)
    assert DeviceConfig(
        epi_thickness_um=1.0, buffer_thickness_um=None, buffer_end_doping=None
    ).buffer_thickness_um is None
    with pytest.raises(ValueError, match="below"):
        DeviceConfig(buffer_thickness_um=10.0, buffer_end_doping=1e16)
    with pytest.raises(ValueError):
        DeviceConfig(buffer_thickness_um=1.0, buffer_end_doping=-1.0)


def test_devsim_donors_equal_the_numpy_twin():
    name = f"mcv_{uuid.uuid4().hex[:6]}"
    info = create_dd_device(device_name=name, **_preset_kwargs())
    try:
        g = lambda n: np.array(
            devsim.get_node_model_values(device=name, region="sic", name=n)
        )
        x, d = g("x"), g("Donors")
    finally:
        devsim.delete_device(device=name)
    assert info["buffer_thickness_cm"] == M["buffer_thickness_cm"]
    twin = donor_profile(
        x, SUB, M["N_D_junction"], M["N_D_bulk"], M["L_transition"], buffer=BUFFER
    ) + 1e19 * (x >= SUB + EPI)
    nz = twin > 0
    assert np.max(np.abs(d[nz] - twin[nz]) / twin[nz]) < 1e-13
    # the transition: epi value at its start, N_end just before the epi end, 1e19 beyond
    x0 = SUB + EPI - M["buffer_thickness_cm"]
    i0, i1 = np.searchsorted(x, x0), np.searchsorted(x, SUB + EPI)
    assert (
        abs(d[i0] / 3.5e13 - 1) < 0.05
        and abs(d[i1 - 1] / M["buffer_end_doping"] - 1) < 0.01
    )
    assert d[i1] > 1e19


# ------------------------------------------------------------ parity engine vs devsim (1D)
BIASES = np.array([0.0, -2.0, -10.0, -40.0, -100.0])


def test_engine_matches_devsim_with_buffer():
    kw = _preset_kwargs()
    name = f"mcvp_{uuid.uuid4().hex[:6]}"
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
        M["N_D_junction"],
        M["N_D_bulk"],
        M["L_transition"],
        nplus=(SUB + EPI, 1e19),
        buffer=BUFFER,
    )
    res = reverse_sweep(pb, BIASES, SUB, EPI, keep_psi=True)
    assert res.converged.all()
    np.testing.assert_allclose(res.W[:, 0], ref["depletion_widths"], rtol=1e-10)
    np.testing.assert_allclose(res.psi[-1, 0], psi_ref, atol=1e-9)
    c = cce_vs_bias(BIASES, device_kwargs=kw)
    G = alpha_generation(x, SUB, 15e-4, 1e18, epi_end=SUB + EPI)
    ce = cce_sweep(pb, BIASES, SUB, G)
    assert ce.converged.all()
    np.testing.assert_allclose(ce.cce[:, 0], c["cce_values"], atol=1e-9)


# ------------------------------------------------------------------- engine boxes
def test_preset_inside_every_engine_domain():
    from etna.core.engine2d import (
        engine2d_dark_unsupported_reason,
        engine2d_unsupported_reason,
    )

    c = measured_cv_config()
    assert sim.validated_domain_reason(c) is None
    kw = device_kwargs(c)
    kw.pop("N_D"), kw.pop("epi_thickness_cm")
    nd = kw.pop("N_D_bulk")
    for hw in (25.0, 50.0):
        for vb in (20.0, 40.0, 60.0):
            assert engine2d_unsupported_reason(hw, vb, 10e-4, nd, 10, **kw) is None
    assert engine2d_dark_unsupported_reason(50.0, 50.0, N_D_bulk=nd, **kw) is None
    # box edges
    assert sim._ENGINE_BOX["N_D_bulk"][0] == 2e13
    assert sim.validated_domain_reason(DeviceConfig(N_D_bulk=2e13)) is None
    assert "N_D_bulk" in sim.validated_domain_reason(DeviceConfig(N_D_bulk=1.9e13))


@pytest.mark.parametrize(
    "over, fragment",
    [
        (dict(buffer_end_doping=5e15), "buffer_end_doping"),
        (dict(buffer_thickness_um=3.0), "buffer_thickness_um"),
        (dict(buffer_thickness_um=0.3), "buffer_thickness_um"),
    ],
)
def test_buffer_outside_box_goes_to_devsim(over, fragment):
    from etna.core.engine2d import engine2d_unsupported_reason

    c = measured_cv_config(**over)
    assert fragment in sim.validated_domain_reason(c)
    kw = device_kwargs(c)
    kw.pop("N_D"), kw.pop("epi_thickness_cm")
    nd = kw.pop("N_D_bulk")
    # the 2D engine names the cm keyword
    assert fragment.replace("_um", "") in engine2d_unsupported_reason(
        25.0, 20.0, 10e-4, nd, 3, **kw
    )


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
    """Plateau (-5..-59 V) rms <= 0.02 um and the three labels; the knee only to 0.3 um rms:
    the exponential epi cannot make the 1.3e13 dip at 4-6 um (LOG.md), so the 13 knee
    markers are NOT within their digitisation sigma (max 3.9 sigma) -- a stated limit of the
    parametric form, not of the implementation (the 19-knot profile would be needed)."""
    cfg = measured_cv_config(area_cm2=0.04)
    V = np.round(np.arange(0.0, -60.001, -0.25), 2)
    r = sim.run_cv(cfg, 0.0, -60.0, len(V))
    assert r.metadata["solver"] == "batched" and r.metadata["validated_domain"]
    W = EPS * cfg.area_cm2 / r.y * 1e4
    at = lambda v: W[np.argmin(np.abs(r.x - v))]
    # the lab fit's own numbers (fits.json fit_buffer_Nnp), rounded record
    np.testing.assert_allclose(
        [at(0), at(-10), at(-30)], [1.6547, 9.5185, 9.7503], atol=2e-4
    )
    # labels 1.68 / 9.52 / 9.73
    assert (
        abs(at(0) - 1.68) < 0.03
        and abs(at(-10) - 9.52) < 0.01
        and abs(at(-30) - 9.73) < 0.025
    )
    plateau, knee = _digitised()
    res_p = np.array([w - at(v) for v, w in plateau])
    assert len(plateau) > 200 and np.sqrt(np.mean(res_p**2)) <= 0.02
    res_k = np.array([w - at(v) for v, w in knee])
    assert len(knee) == 13 and np.sqrt(np.mean(res_k**2)) < 0.3
    assert (
        abs(at(-2) - 7.18) < 0.3 and abs(at(-5) - 9.27) < 0.03
    )  # default: 2.99 / 5.93
    # full depletion of the low-doped epi (edge within L_D of the transition start)
    assert r.metadata["reach_through_V"] == pytest.approx(3.1, abs=0.051)


def test_analytic_vfd_profile_integral_matches_closed_form_without_buffer():
    from etna.export.labels import analytic_vfd
    from etna.core.analytical import (
        full_depletion_voltage_graded,
        full_depletion_voltage_profile,
    )

    c = nplus_calibrated_config()
    y = np.linspace(0.0, c.epi_thickness_um * 1e-4, 20001)
    nd = donor_profile(y, 0.0, c.N_D_junction, c.N_D_bulk, c.L_transition_um * 1e-4)
    closed = full_depletion_voltage_graded(
        c.epi_thickness_um * 1e-4,
        c.N_D_bulk,
        c.N_D_junction,
        c.L_transition_um * 1e-4,
        0.0,
    )
    assert full_depletion_voltage_profile(y, nd, 0.0) == pytest.approx(closed, rel=1e-6)
    assert analytic_vfd(c) == pytest.approx(9.93, abs=0.02)  # the previous default
    # with the transition V_fd is the depletion of the low-doped epi up to the transition
    # start (8.86 um): 2.74 V, next to the engine's 3.1 V reach-through; integrating through
    # the transition to 10 um would give ~130 V (the C-V creep, not full depletion)
    v_pre = analytic_vfd(measured_cv_config())
    assert v_pre == pytest.approx(2.736, abs=0.01)
    v_nob = analytic_vfd(
        measured_cv_config(buffer_thickness_um=None, buffer_end_doping=None)
    )
    assert 0 < v_pre < v_nob < 5.0


# ------------------------------------------------------------------ exports and the app
def test_spec_sheet_and_cross_section_label_the_preset():
    from etna.export.cross_section import cross_section_svg, stack_layers
    from etna.export.labels import MEASURED_CV_NOTE, epi_donor_profile
    from etna.export.spec_sheet import build_spec, spec_to_markdown

    cfg = measured_cv_config(area_cm2=0.04)
    spec = build_spec(cfg, None)
    assert spec["measured_cv_preset"] and spec["reference_device"]
    assert "measured-C-V profile" in spec["model"] and "default stack" in spec["model"]
    assert spec["cv"]["calibration"] == MEASURED_CV_NOTE
    cal = spec["cv"]["calibration"]
    # the statements the default's export must carry
    for s in ("digitised Photons Fig. 4b", "plateau rms 0.013 um", "knee rms 0.24 um",
              "below the documented 0.5-1e14", "not by SIMS"):
        assert s in cal, s
    assert "fit input" in spec["statements"][1].lower()
    names = [L["name"] for L in stack_layers(cfg)]
    assert names == ["p+ anode", "n- epi", "n-/n+ transition", "n+ substrate"]
    assert "transition" in cross_section_svg(cfg) and "transition" in spec_to_markdown(
        spec
    )
    assert float(epi_donor_profile(cfg, 10.0)) == pytest.approx(
        M["buffer_end_doping"], rel=1e-3
    )
    assert float(epi_donor_profile(cfg, 5.0)) == pytest.approx(M["N_D_bulk"], rel=1e-2)
    # the previous default is labelled as such, with its 3-point calibration note
    d = build_spec(nplus_calibrated_config(), None)
    assert d["model"].startswith("previous default") and not d["reference_device"]
    assert d["nplus_calibrated_preset"] and "0 V and -10 V" in d["cv"]["calibration"]
    assert d["epi_doping"]["buffer_thickness_um"] is None


def test_sidebar_values_with_the_preset_assemble_the_preset_config():
    from app.components.device_sidebar import assemble_config

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
    # without the transition keys the sidebar assembles the default transition
    assert assemble_config(values).buffer_thickness_um == 1.136634
    values.update(etna.MEASURED_CV_PRESET)
    assert assemble_config(values) == measured_cv_config()  # 0.04 cm2 = the default area
    values.update(etna.NPLUS_CALIBRATED_PRESET)
    assert assemble_config(values) == nplus_calibrated_config()


def test_run_all_sweep_passes_the_buffer_to_the_engine():
    rep = etna.run_all(measured_cv_config(area_cm2=0.04), sections=("sweep",))
    s = rep.get("sweep")
    assert s.status == "ok" and s.solver == "engine2d", (s.status, s.reason)
    assert "transition" in s.summary
    df = s.result.metadata["table"]
    assert df["V_fd"].iloc[0] == pytest.approx(
        __import__("etna.export.labels", fromlist=["analytic_vfd"]).analytic_vfd(
            measured_cv_config()
        ),
        rel=1e-3,
    )


# ---------------------------------------------------------------------------- slow
@pytest.mark.slow
def test_engine2d_matches_devsim_with_buffer():
    from etna.core.charge_collection_2d import cce_lateral_scan, create_2d_dd_device
    from etna.core.engine2d import lateral_cce_scan_2d

    kw = _preset_kwargs()
    kw.pop("N_D"), kw.pop("doping_profile")
    e = lateral_cce_scan_2d(25.0, 20.0, n_lateral_points=3, **kw)
    name = f"mcv2d_{uuid.uuid4().hex[:8]}"
    info = create_2d_dd_device(half_width_um=25.0, V_bias=20.0, device_name=name, **kw)
    try:
        d = cce_lateral_scan(info, n_points=3)
    finally:
        devsim.delete_device(device=name)
    np.testing.assert_allclose(e["cce_values"], d["cce_values"], rtol=0, atol=1e-10)


@pytest.mark.slow
def test_nt_reproduces_the_18pA_bound_on_the_preset():
    from etna.core.dark_current import (
        AREA_DOSIMETER_CM2,
        create_dark_current_device,
        dark_current_sweep,
        setup_tat_model,
    )

    info = create_dark_current_device(
        N_t=M["N_t"],
        E_t=M["E_t"],
        device_name=f"mcvnt_{uuid.uuid4().hex[:6]}",
        **_preset_kwargs(),
    )
    try:
        setup_tat_model(info, E_t=M["E_t"], N_t=M["N_t"])
        r = dark_current_sweep(info, np.array([-60.0]), area=AREA_DOSIMETER_CM2)
    finally:
        devsim.delete_device(device=info["device_name"])
    assert abs(r["I_total"][0]) == pytest.approx(18e-12, rel=1e-5)

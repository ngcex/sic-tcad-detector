"""Pins the n+ calibration record (etna.core.device.NPLUS_CALIBRATION): the default model
2026-10-17 .. 2026-10-18, now the named preset etna.nplus_calibrated_config() (the default
is the measured-C-V profile, DEFAULT_STACK). Abrupt epi/n+ step: no graded transition.

C-V (since 2026-10-17): the targets are W = eps*A/C of a measured C, so the model
observable is eps/C of the small-signal C = dQ/dV (reverse_sweep's C), not the depletion
edge. 0 V and -10 V are solved exactly at the held L_transition (calibration inputs); the
-30 V target (9.73 um) is out of fit (+2.56 %), pinned too. N_t: I(300 K, -60 V, 4 mm2) =
18 pA. See lab/solver2d/bench/nplus_calibrate.py (cal-c).
"""

import uuid

import devsim
import numpy as np
import pytest
import scipy.optimize

from etna.core.batched import make_batch, node_positions, reverse_sweep
from etna.core.batched.material import EPS
from etna.core.device import NPLUS_CALIBRATION as CAL
from etna.core.radiation_damage import RadiationDamageParams
from etna.core.sic_material import SiC4H_Parameters

EPI = 10e-4
V_T = np.array([0.0, -10.0, -30.0])
W_T = np.array([1.7e-4, 9.5e-4, 9.73e-4])


def _W(NJ, NB, L):
    sub, t_np = CAL["substrate_thickness_cm"], CAL["nplus_thickness_cm"]
    x = node_positions(EPI, sub, t_np)
    pb = make_batch(x, sub, 1e19, NJ, NB, L, 300.0, nplus=(sub + EPI, CAL["N_D_nplus"]))
    r = reverse_sweep(pb, V_T, sub, EPI)
    assert r.converged.all()
    return EPS / r.C[:, 0]  # eps*A/C, what the measured W targets are


def test_trap_level_reconciled_in_the_record():
    """The record uses the radiation-damage Z1/2 level (Burin et al. Table I), applied to
    SiC4H_Parameters with the default flip (0.65 -> 0.67). On the no-n+ model it made the
    dark sweep stop at -89.5 V, so the legacy preset keeps 0.65; with n+ the default reaches
    -100 V (tests/test_default_stack.py)."""
    assert CAL["E_t"] == RadiationDamageParams().E_Z12 == 0.67
    assert SiC4H_Parameters().E_t == 0.67


def test_record_solves_the_0V_and_10V_points_exactly():
    W = _W(CAL["N_D_junction"], CAL["N_D_bulk"], CAL["L_transition"])
    np.testing.assert_allclose(W[:2], W_T[:2], rtol=2e-6)  # record rounded to 7 digits
    # -30 V: out of fit, +2.56 % (buffer / thinner epi / fringe: not separable, LOG.md)
    assert (W[2] - W_T[2]) / W_T[2] == pytest.approx(0.0256, abs=5e-4)
    assert (
        0.5e14 <= CAL["N_D_bulk"] <= 1e14
    )  # paper: epi 0.5-1e14 cm^-3 (independent check)


def test_two_point_solve_reproduces_the_record():
    L = CAL["L_transition"]

    def res(z):
        return ((_W(np.exp(z[0]), np.exp(z[1]), L) - W_T) / W_T)[:2]

    r = scipy.optimize.least_squares(
        res, np.log([2.974271e15, 6.335936e13]), xtol=1e-15, ftol=1e-15, gtol=1e-15
    )
    np.testing.assert_allclose(
        np.exp(r.x), [CAL["N_D_junction"], CAL["N_D_bulk"]], rtol=1e-6
    )


def test_three_point_system_has_no_exact_solution():
    """Least squares on all three targets: a residual of ~3 % remains (why the record
    holds L and solves only 0 / -10 V). On the old edge observable N_D_bulk also ran to
    its 1e12 bound; on eps/C it does not (lab/microdosimeter/notes/buffer_w.out)."""
    lo, hi = np.log([1e14, 1e12, 0.05e-4]), np.log([1e17, 1e15, 5e-4])
    r = scipy.optimize.least_squares(
        lambda z: (_W(*np.exp(z)) - W_T) / W_T,
        np.log([2.9e15, 8.5e13, 1e-4]),
        bounds=(lo, hi),
    )
    assert np.abs(r.fun).max() > 0.01


@pytest.mark.slow
def test_nt_reproduces_the_18pA_bound():
    from etna.core.dark_current import (
        AREA_DOSIMETER_CM2,
        create_dark_current_device,
        dark_current_sweep,
    )

    kw = {
        k: CAL[k]
        for k in (
            "substrate_thickness_cm",
            "nplus_thickness_cm",
            "N_D_nplus",
            "N_D_junction",
            "N_D_bulk",
            "L_transition",
            # None: the record has no graded transition (the builders' default since
            # 2026-10-18 would add the measured-C-V one)
            "buffer_thickness_cm",
            "buffer_end_doping",
        )
    }
    info = create_dark_current_device(
        N_t=CAL["N_t"],
        device_name=f"cal_{uuid.uuid4().hex[:6]}",
        doping_profile="graded",
        **kw,
    )
    try:
        from etna.core.dark_current import setup_tat_model

        setup_tat_model(info, E_t=CAL["E_t"], N_t=CAL["N_t"])  # the record's trap level
        r = dark_current_sweep(info, np.array([-60.0]), area=AREA_DOSIMETER_CM2)
    finally:
        devsim.delete_device(device=info["device_name"])
    assert abs(r["I_total"][0]) == pytest.approx(18e-12, rel=1e-5)

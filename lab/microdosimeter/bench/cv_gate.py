"""GATE: independent checks of the small-signal W_C = eps / (dQ/dV) of bench/buffer_w.py.

buffer_w.py derived W_C from the displacement current through a carrier-free plane
x0 = 0.5 um below the junction (C = eps dE(x0)/dV). Here, on the production
DEFAULT_STACK, two independent routes:

(a) batched engine, charge derivatives (no plane): in a reverse-biased diode the
    electrons that leave/enter the device all pass the cathode and the holes the anode
    (generation/recombination negligible in the quasi-static limit), so the terminal
    charge is Q = q * sum_i v_i n_i = q * sum_i v_i p_i + const (box-integration volumes
    v_i, the discretisation of the Poisson solver). Central differences with
    DV = 0.05 / 0.1 / 0.25 V. Also: the space charge of the n side from the plane x0
    (q sum_{x>=x0} v (D - n + p)), and from the metallurgical junction (to show why
    buffer_w.py avoided it).
(b) devsim small-signal AC: the DD device of build_device(DEFAULT config), cathode
    re-attached to a circuit voltage source (acreal = 1), devsim.solve(type="ac"),
    C = Im(I_V1) / (2 pi f) at f = 1 Hz and 1 kHz.

usage: uv run python lab/microdosimeter/bench/cv_gate.py
"""

import math

import numpy as np

from etna.core.batched import make_batch, node_positions
from etna.core.batched.extract import depletion_width
from etna.core.batched.material import EPS, Q
from etna.core.batched.sweeps import _Ramp
from etna.core.device import DEFAULT_STACK as S

V_T = np.array([0.0, -10.0, -30.0])
W_LAB_C = np.array([1.708, 9.764, 9.981])  # um, lab W_C (notes/buffer_w.out)
EPI, SUB, NP = 10e-4, S["substrate_thickness_cm"], S["nplus_thickness_cm"]
X0 = 0.5e-4


def engine_states(volts, tol=1e-12):
    x = node_positions(EPI, SUB, NP)
    pb = make_batch(
        x,
        SUB,
        1e19,
        S["N_D_junction"],
        S["N_D_bulk"],
        S["L_transition"],
        300.0,
        nplus=(SUB + EPI, S["N_D_nplus"]),
        buffer=None,
    )
    ramp = _Ramp(pb, tol)
    out = {}
    # walk from the most forward point to the deepest reverse point
    for V in sorted(volts, reverse=True):
        res = ramp.to(-V)
        assert res.converged.all(), V
        out[V] = res
    return pb, out


def main():
    dvs = (0.05, 0.1, 0.25)
    volts = sorted(
        {round(v + s * d, 9) for v in V_T for d in dvs for s in (-1, 0, 1)},
        reverse=True,
    )
    pb, st = engine_states(volts)
    x = pb.x
    h = np.diff(x)
    vol = np.empty_like(x)
    vol[1:-1] = 0.5 * (h[1:] + h[:-1])
    vol[0], vol[-1] = 0.5 * h[0], 0.5 * h[-1]
    xm = 0.5 * (x[1:] + x[:-1])
    D = pb.netdop[0]

    def Qs(res):
        n, p, psi = res.n[0], res.p[0], res.psi[0]
        E = -np.diff(psi) / h
        sel0 = x >= SUB + X0
        selj = x >= SUB
        return dict(
            electrons=Q * np.sum(vol * n),
            holes=Q * np.sum(vol * p),
            nside_x0=Q * np.sum((vol * (D - n + p))[sel0]),
            nside_junction=Q * np.sum((vol * (D - n + p))[selj]),
            field_x0=EPS * np.interp(SUB + X0, xm, E),
        )

    print("== (a) batched engine, DEFAULT_STACK, W = eps / (dQ/dV) in um ==")
    print(f"lab W_C (buffer_w.py, plane x0 = 0.5 um, DV 0.25): {W_LAB_C}")
    for V in V_T:
        We = depletion_width(x, st[V].n, pb.donors, SUB, EPI + NP)[0]
        print(f"V = {V:+.0f} V   W_edge = {We*1e4:.4f}")
        for key in ("electrons", "holes", "nside_x0", "nside_junction", "field_x0"):
            row = []
            for dv in dvs:
                lo, hi = round(V - dv, 9), round(V + dv, 9)
                dQ = Qs(st[lo])[key] - Qs(st[hi])[key]
                row.append(EPS * 2 * dv / abs(dQ) * 1e4)
            print(
                f"   {key:<15}"
                + "".join(f"  DV {dv:4.2f}: {w:8.4f}" for dv, w in zip(dvs, row))
            )
    return pb


def devsim_ac(volts=(0.0, -10.0, -30.0), freqs=(1.0, 1e3), cfg=None):
    """devsim AC C (F/cm^2) per bias, {V: C at freqs[0]}. Run it in its own process and
    end the process afterwards: devsim 2.10 cannot remove a circuit again
    (delete_circuit / reset_devsim corrupt the next DC solve)."""
    import devsim
    import devsim.python_packages.simple_physics as sp

    from etna.api.device import DeviceConfig, build_device
    from etna.core.devsim_reset import reset_devsim_fully
    from etna.core.poisson import extract_depletion_width_numerical

    def dc():
        # the circuit row converges to ~1e-11 relative only: 1e-10, else the production 1e-8
        try:
            devsim.solve(
                type="dc", absolute_error=1e10, relative_error=1e-10, maximum_iterations=60
            )
        except devsim.error:
            devsim.solve(
                type="dc", absolute_error=1e10, relative_error=1e-8, maximum_iterations=60
            )

    reset_devsim_fully()
    cfg = DeviceConfig() if cfg is None else cfg
    info = build_device(cfg)
    dev, reg = info["device_name"], info["region_name"]
    bias = sp.GetContactBiasName("cathode")
    devsim.circuit_element(name="V1", n1=bias, n2=0, value=0.0, acreal=1.0, acimag=0.0)
    sp.CreateSiliconPotentialOnlyContact(dev, reg, "cathode", True)
    sp.CreateSiliconDriftDiffusionAtContact(dev, reg, "cathode", True)
    dc()
    print("== (b) devsim small-signal AC (cathode on V1, acreal = 1) ==")
    Vc_now = 0.0
    out = {}
    for V in volts:
        target = -V
        steps = max(1, int(math.ceil(abs(target - Vc_now) / 0.5)))
        for v in np.linspace(Vc_now, target, steps + 1)[1:]:
            # the contact model reads the device parameter (it shadows the circuit
            # node in the DC residual); the circuit node carries the AC excitation
            devsim.circuit_alter(name="V1", value=float(v))
            devsim.set_parameter(device=dev, name=bias, value=float(v))
            dc()
        Vc_now = target
        row = []
        for f in freqs:
            devsim.solve(type="ac", frequency=f)
            Iim = devsim.get_circuit_node_value(node="V1.I", solution="ssac_imag")
            Ire = devsim.get_circuit_node_value(node="V1.I", solution="ssac_real")
            C = abs(Iim) / (2 * math.pi * f)  # F/cm^2 (1D)
            row.append((f, EPS / C * 1e4, Ire))
            out.setdefault(float(V), C)
        We = extract_depletion_width_numerical(info)
        print(
            f"V = {V:+.0f} V  W_edge {We*1e4:.4f}  "
            + "  ".join(
                f"f {f:g} Hz: W_C {w:.7f} um (Re I {re:.2e})" for f, w, re in row
            )
        )
    return out


if __name__ == "__main__":
    main()
    devsim_ac()

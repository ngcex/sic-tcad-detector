"""1D batched engine: mesh convergence of CCE (run_cce's alpha generation), small-signal
C and the depletion edge W on the production mesh (devsim's build_diode_mesh = level 0)
under (a) uniform bisection of every interval and (b) bisection of the two contact
cells only. Richardson extrapolation from the last three uniform levels with the
observed order.

Two CCE normalisations are reported: the engine's (and devsim 1D's) J / (q trapz(G)),
which counts the generation of the Dirichlet contact nodes' half cells, and the
"injected" one of the 2D path (compute_cce_2d since 2026-10-18), J / (q sum_i G_i V_i)
over the non-contact nodes only. They differ only where G is nonzero on a contact node.

Usage: uv run python lab/mesh_contact/bench/default_1d.py [legacy] [uniform]
  default: DeviceConfig() (n+ substrate, graded transition)
  legacy:  DeviceConfig(nplus_thickness_um=None, buffer_*=None): cathode on the n- epi
  uniform: uniform generation over the epi (FLASH-like) instead of the alpha profile
"""

import sys
import time

import numpy as np

np.set_printoptions(floatmode="maxprec", formatter={"float": lambda v: f"{v:.3e}"})

from etna import DeviceConfig
from etna.api.device import device_kwargs
from etna.api.models import refine_mesh
from etna.api.simulation import _engine_batch
from etna.core.batched.mesh import node_positions
from etna.core.batched.sweeps import alpha_generation, cce_sweep, reverse_sweep
from etna.core.generation_profiles import ALPHA_RANGE_SIC_CM

V = np.array([-1.0, -2.0, -3.0, -5.0, -10.0, -30.0, -50.0])


def contact_refine(x, times):
    """Bisect ``times`` times the half of each contact cell touching the contact."""
    x = list(np.asarray(x, float))
    out = set(x)
    for c, o in ((x[0], x[1]), (x[-1], x[-2])):
        for _ in range(times):
            o = 0.5 * (c + o)
            out.add(o)
    return np.array(sorted(out))


UNIFORM = "uniform" in sys.argv


def solve(cfg, x, voltages=V, injected=False):
    kw = device_kwargs(cfg)
    pb, x = _engine_batch([kw], x)
    sub = kw["substrate_thickness_cm"]
    epi = kw["epi_thickness_cm"]
    epi_end = sub + epi if kw.get("nplus_thickness_cm") else None
    if UNIFORM:
        G = np.where((x >= sub) & (x <= (epi_end or x[-1]) + 1e-12), 1e18, 0.0)
    else:
        G = alpha_generation(x, sub, ALPHA_RANGE_SIC_CM, 1e18, epi_end=epi_end)
    r = cce_sweep(pb, voltages, sub, G, T=[kw["T"]], N_A=[kw["N_A"]])
    d = reverse_sweep(pb, voltages, sub, epi)
    assert r.converged.all() and d.converged.all()
    if injected:
        vol = np.zeros_like(x)
        vol[1:] += 0.5 * np.diff(x)
        vol[:-1] += 0.5 * np.diff(x)
        # same collected current, denominator over the non-contact nodes only
        return r.J_collected[:, 0] / (r.J_generated * np.sum((G * vol)[1:-1])
                                      / np.trapezoid(G, x))  # fmt: skip
    return r.cce[:, 0], d.C[:, 0], d.W[:, 0]


def richardson(a0, a1, a2):
    """Observed order p, extrapolated value and |ext - a2| from three levels h, h/2, h/4."""
    d1, d2 = a1 - a0, a2 - a1
    with np.errstate(divide="ignore", invalid="ignore"):
        q = d1 / d2
        p = np.log2(np.abs(q))
        ext = a2 + d2 / (np.abs(q) - 1)
    return p, ext, np.abs(ext - a2)


def legacy_config():
    return DeviceConfig(
        nplus_thickness_um=None, buffer_thickness_um=None, buffer_end_doping=None
    )


if __name__ == "__main__":
    legacy = "legacy" in sys.argv
    cfg = legacy_config() if legacy else DeviceConfig()
    kw = device_kwargs(cfg)
    x0 = node_positions(
        kw["epi_thickness_cm"],
        kw["substrate_thickness_cm"],
        kw.get("nplus_thickness_cm"),
    )
    print("stack:", "LEGACY (cathode on epi)" if legacy else "DEFAULT", " N0 =", len(x0),
          " anode cell %.1f nm, cathode cell %.1f nm" % ((x0[1] - x0[0]) * 1e7,
                                                         (x0[-1] - x0[-2]) * 1e7))  # fmt: skip
    print("V:", V)
    for mode in ("uniform", "contact"):
        res, res_inj = [], []
        for lev in range(5 if mode == "uniform" else 8):
            if lev == 0:
                x = x0
            elif mode == "uniform":
                x = refine_mesh(x0, lev)
            else:
                x = contact_refine(x0, lev)
            t = time.perf_counter()
            c, C, W = solve(cfg, np.asarray(x))
            if mode == "contact":
                ci = solve(cfg, np.asarray(x), injected=True)
                res_inj.append(ci)
                print(f"   injected-normalised CCE - L0: "
                      f"{np.array2string(ci - res_inj[0], precision=2)}  "
                      f"(trapz - injected: {np.array2string(c - ci, precision=2)})")  # fmt: skip
            res.append((c, C, W))
            print(f"{mode:8s} L{lev} N={len(x):6d} ({time.perf_counter() - t:.1f}s) "
                  f"CCE-CCE0={np.array2string(c - res[0][0], precision=2)} "
                  f"dC/C0={np.array2string(C / res[0][1] - 1, precision=2)}", flush=True)  # fmt: skip
        if mode == "uniform":
            for k, name in enumerate(("CCE", "C", "W")):
                a = [r[k] for r in res]
                p, ext, err = richardson(a[-3], a[-2], a[-1])
                e0 = (a[0] - ext) if name == "CCE" else (a[0] / ext - 1)
                ee = err if name == "CCE" else err / np.abs(ext)
                print(f"  {name}: level0 {np.array2string(a[0], precision=8)}")
                print(f"      extrap {np.array2string(ext, precision=8)}")
                print(f"      order  {np.array2string(p, precision=2)}")
                print(f"      level0 error ({'abs' if name == 'CCE' else 'rel'}) "
                      f"{np.array2string(e0, precision=2)}  |ext - L4| "
                      f"{np.array2string(ee, precision=1)}")  # fmt: skip

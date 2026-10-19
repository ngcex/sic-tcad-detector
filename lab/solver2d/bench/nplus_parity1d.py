"""Parity of the batched 1D engine vs devsim WITH an n+ layer (W/C, potential, CCE).

Same mesh (node_positions with nplus_thickness_cm == build_diode_mesh), same doping,
devsim path = cv_sweep / cce_vs_bias with the n+ device keywords.

usage: uv run python lab/solver2d/bench/nplus_parity1d.py
"""

from __future__ import annotations

import uuid

import devsim
import numpy as np

from etna.core.batched import (
    alpha_generation,
    cce_sweep,
    make_batch,
    node_positions,
    reverse_sweep,
)
from etna.core.charge_collection import cce_vs_bias
from etna.core.cv_analysis import cv_sweep
from etna.core.drift_diffusion import create_dd_device

GRADED = dict(
    doping_profile="graded", N_D_junction=2.9e15, N_D_bulk=8.5e13, L_transition=1e-4
)
CV_BIAS = np.array([0.0, -5.0, -10.0, -20.0, -40.0, -60.0, -80.0, -100.0])
CCE_BIAS = np.array([0.0, -5.0, -10.0, -30.0, -60.0, -100.0])

CASES = [
    ("default np 1um 1e19", dict(**GRADED), (1e-4, 1e19)),
    ("default np 1um 1e18", dict(**GRADED), (1e-4, 1e18)),
    ("default np 5um 1e19", dict(**GRADED), (5e-4, 1e19)),
    ("p+ 0.3um np 1um", dict(**GRADED, substrate_thickness_cm=0.3e-4), (1e-4, 1e19)),
    (
        "hot thin N_A 5e18",
        dict(
            **GRADED,
            T=350.0,
            N_A=5e18,
            epi_thickness_cm=8e-4,
            substrate_thickness_cm=2e-4,
        ),
        (1e-4, 1e19),
    ),
    ("cold 250K", dict(**GRADED, T=250.0), (1e-4, 1e19)),
    ("400K", dict(**GRADED, T=400.0), (1e-4, 1e19)),
    (
        "epi 5um bulk 5e13",
        dict(
            doping_profile="graded",
            N_D_junction=1e15,
            N_D_bulk=5e13,
            L_transition=0.3e-4,
            epi_thickness_cm=5e-4,
        ),
        (1e-4, 1e19),
    ),
    (
        "epi 30um bulk 2e15",
        dict(
            doping_profile="graded",
            N_D_junction=1e16,
            N_D_bulk=2e15,
            L_transition=3e-4,
            epi_thickness_cm=30e-4,
            substrate_thickness_cm=3e-4,
        ),
        (1e-4, 1e19),
    ),
    ("uniform 1e15", dict(doping_profile="uniform", N_D=1e15), (1e-4, 1e19)),
]


def batch_for(kw, np_):
    epi = kw.get("epi_thickness_cm", 10e-4)
    sub = kw.get("substrate_thickness_cm", 1e-4)
    if kw["doping_profile"] == "uniform":
        nj = nb = kw["N_D"]
        L = 1e-4
    else:
        nj, nb, L = kw["N_D_junction"], kw["N_D_bulk"], kw["L_transition"]
    x = node_positions(epi, sub, np_[0])
    pb = make_batch(
        x,
        sub,
        kw.get("N_A", 1e19),
        nj,
        nb,
        L,
        kw.get("T", 300.0),
        nplus=(sub + epi, np_[1]),
        buffer=None,
    )
    return x, pb, epi, sub


def main():
    worst = dict(W=0.0, psi=0.0, cce=0.0)
    for label, kw, np_ in CASES:
        dkw = dict(kw, nplus_thickness_cm=np_[0], N_D_nplus=np_[1])
        name = f"par_{uuid.uuid4().hex[:6]}"
        info = create_dd_device(device_name=name, **dkw)
        ref = cv_sweep(info, CV_BIAS)
        psi_ref = np.array(
            devsim.get_node_model_values(device=name, region="sic", name="Potential")
        )
        x_ref = np.array(
            devsim.get_node_model_values(device=name, region="sic", name="x")
        )
        devsim.delete_device(device=name)
        x, pb, epi, sub = batch_for(kw, np_)
        assert np.array_equal(x, x_ref), label
        res = reverse_sweep(pb, CV_BIAS, sub, epi, keep_psi=True)
        nref = len(ref["depletion_widths"])
        dW = np.max(
            np.abs(res.W[:nref, 0] - ref["depletion_widths"]) / ref["depletion_widths"]
        )
        dpsi = np.max(np.abs(res.psi[nref - 1, 0] - psi_ref))
        c = cce_vs_bias(CCE_BIAS, epi_thickness_cm=epi, device_kwargs=dkw)
        G = alpha_generation(x, sub, 15e-4, 1e18, epi_end=sub + epi)
        ce = cce_sweep(
            pb,
            CCE_BIAS,
            sub,
            G,
            T=kw.get("T", 300.0),
            N_A=kw.get("N_A", 1e19),
        )
        nc = len(c["cce_values"])
        dc = np.max(np.abs(ce.cce[:nc, 0] - c["cce_values"]))
        raw_max = np.max(ce.J_collected[:, 0] / ce.J_generated)
        worst["W"], worst["psi"], worst["cce"] = (
            max(worst["W"], dW),
            max(worst["psi"], dpsi),
            max(worst["cce"], dc),
        )
        print(
            f"RES {label:22s} cv_pts {nref}/{len(CV_BIAS)} trunc {ref['truncated']} conv {res.converged.all()} "
            f"max rel dW {dW:.2e} max |dpsi| {dpsi:.2e} V | cce_pts {nc}/{len(CCE_BIAS)} conv {ce.converged.all()} "
            f"max |dCCE| {dc:.2e} max raw CCE {raw_max:.6f} | CCE(-100) {ce.cce[-1, 0]:.6f}",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()

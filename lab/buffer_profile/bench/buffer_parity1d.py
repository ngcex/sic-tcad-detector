"""Parity of the batched 1D engine vs devsim WITH the graded n-/n+ transition (buffer),
and at the low N_D_bulk of the measured-C-V preset (box extension 5e13 -> 2e13).

Same mesh (node_positions == build_diode_mesh), same doping expression (devsim string vs
its numpy twin etna.core.device.donor_profile), devsim path = cv_sweep / cce_vs_bias with
the buffer keywords. Reports max rel dW, max |dpsi| (V) and max |dCCE| per case.

usage: uv run python lab/buffer_profile/bench/buffer_parity1d.py
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
from etna.core.device import MEASURED_CV_STACK as M
from etna.core.device import mobility_doping
from etna.core.drift_diffusion import create_dd_device

CV_BIAS = np.array([0.0, -2.0, -5.0, -10.0, -20.0, -40.0, -60.0, -80.0, -100.0])
CCE_BIAS = np.array([0.0, -2.0, -5.0, -10.0, -30.0, -60.0, -100.0])

PRESET = dict(
    doping_profile="graded",
    N_D_junction=M["N_D_junction"],
    N_D_bulk=M["N_D_bulk"],
    L_transition=M["L_transition"],
    substrate_thickness_cm=M["substrate_thickness_cm"],
    nplus_thickness_cm=M["nplus_thickness_cm"],
    N_D_nplus=M["N_D_nplus"],
)
BUF = dict(
    buffer_thickness_cm=M["buffer_thickness_cm"],
    buffer_end_doping=M["buffer_end_doping"],
)
DEF = dict(
    doping_profile="graded",
    N_D_junction=2.958237e15,
    N_D_bulk=8.252732e13,
    L_transition=0.987e-4,
    substrate_thickness_cm=0.3e-4,
    nplus_thickness_cm=1e-4,
    N_D_nplus=1e19,
)

CASES = [
    ("measured-CV preset (buffer)", dict(PRESET, **BUF)),
    ("preset doping, no buffer (N_D_bulk 3.4e13)", dict(PRESET)),
    ("default stack, N_D_bulk 2e13 (box corner)", dict(DEF, N_D_bulk=2e13)),
    (
        "default stack + buffer 0.5um 1e16",
        dict(DEF, buffer_thickness_cm=0.5e-4, buffer_end_doping=1e16),
    ),
    (
        "default stack + buffer 2.5um 1e19",
        dict(DEF, buffer_thickness_cm=2.5e-4, buffer_end_doping=1e19),
    ),
    (
        "bulk 2e13 + buffer 2.5um 1e16",
        dict(DEF, N_D_bulk=2e13, buffer_thickness_cm=2.5e-4, buffer_end_doping=1e16),
    ),
    ("preset, 250 K", dict(PRESET, **BUF, T=250.0)),
    ("preset, 400 K", dict(PRESET, **BUF, T=400.0)),
    (
        "preset, epi 5um buffer 0.5um",
        dict(
            PRESET,
            epi_thickness_cm=5e-4,
            buffer_thickness_cm=0.5e-4,
            buffer_end_doping=M["buffer_end_doping"],
        ),
    ),
    (
        "preset, epi 30um buffer 2.5um n+ 5um",
        dict(
            PRESET,
            epi_thickness_cm=30e-4,
            nplus_thickness_cm=5e-4,
            buffer_thickness_cm=2.5e-4,
            buffer_end_doping=1e18,
        ),
    ),
    (
        "uniform 1e14 + buffer 1um 1e17",
        dict(
            doping_profile="uniform",
            N_D=1e14,
            substrate_thickness_cm=0.3e-4,
            nplus_thickness_cm=1e-4,
            N_D_nplus=1e19,
            buffer_thickness_cm=1e-4,
            buffer_end_doping=1e17,
        ),
    ),
]


def batch_for(kw):
    epi = kw.get("epi_thickness_cm", 10e-4)
    sub = kw["substrate_thickness_cm"]
    t_np = kw["nplus_thickness_cm"]
    if kw["doping_profile"] == "uniform":
        nj = nb = kw["N_D"]
        L = 1e-4
    else:
        nj, nb, L = kw["N_D_junction"], kw["N_D_bulk"], kw["L_transition"]
    x = node_positions(epi, sub, t_np)
    buffer = None
    if kw.get("buffer_thickness_cm") is not None:
        buffer = (sub + epi, kw["buffer_thickness_cm"], kw["buffer_end_doping"])
    pb = make_batch(
        x,
        sub,
        kw.get("N_A", 1e19),
        nj,
        nb,
        L,
        kw.get("T", 300.0),
        nplus=(sub + epi, kw["N_D_nplus"]),
        buffer=buffer,
    )
    return x, pb, epi, sub


def main():
    worst = dict(W=0.0, psi=0.0, cce=0.0, donors=0.0)
    for label, kw in CASES:
        name = f"bpar_{uuid.uuid4().hex[:6]}"
        info = create_dd_device(device_name=name, **kw)
        ref = cv_sweep(info, CV_BIAS)
        g = lambda n: np.array(
            devsim.get_node_model_values(device=name, region="sic", name=n)
        )
        psi_ref, x_ref, don_ref = g("Potential"), g("x"), g("Donors")
        devsim.delete_device(device=name)
        x, pb, epi, sub = batch_for(kw)
        assert np.array_equal(x, x_ref), label
        nz = don_ref > 0
        ddon = np.max(np.abs(pb.donors[0, nz] - don_ref[nz]) / don_ref[nz])
        res = reverse_sweep(pb, CV_BIAS, sub, epi, keep_psi=True)
        nref = len(ref["depletion_widths"])
        dW = np.max(
            np.abs(res.W[:nref, 0] - ref["depletion_widths"]) / ref["depletion_widths"]
        )
        dpsi = np.max(np.abs(res.psi[nref - 1, 0] - psi_ref))
        c = cce_vs_bias(CCE_BIAS, epi_thickness_cm=epi, device_kwargs=kw)
        G = alpha_generation(x, sub, 15e-4, 1e18, epi_end=sub + epi)
        ce = cce_sweep(
            pb,
            CCE_BIAS,
            sub,
            G,
            T=kw.get("T", 300.0),
            N_A=kw.get("N_A", 1e19),
            N_D=mobility_doping(
                kw.get("N_D"), kw["doping_profile"], kw.get("N_D_bulk")
            ),
        )
        nc = len(c["cce_values"])
        dc = np.max(np.abs(ce.cce[:nc, 0] - c["cce_values"]))
        for k, v in (("W", dW), ("psi", dpsi), ("cce", dc), ("donors", ddon)):
            worst[k] = max(worst[k], v)
        print(
            f"RES {label:44s} cv {nref}/{len(CV_BIAS)} trunc {ref['truncated']} conv {res.converged.all()} "
            f"dDonors {ddon:.1e} dW {dW:.2e} dpsi {dpsi:.2e} V | cce {nc}/{len(CCE_BIAS)} conv {ce.converged.all()} "
            f"dCCE {dc:.2e} CCE(-100) {ce.cce[-1, 0]:.6f}",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()

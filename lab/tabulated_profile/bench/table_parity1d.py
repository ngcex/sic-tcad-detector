"""Parity of the batched 1D engine vs devsim with a TABULATED epi donor profile.

Same mesh (node_positions == build_diode_mesh); the donors are the same numpy function in
both paths (etna.core.device.tabulated_donor_profile: the devsim builder sets its values as
a node solution, the engine evaluates it), so the first check is bit-equality of the doping,
then W (edge), psi and CCE through cv_sweep / cce_vs_bias vs reverse_sweep / cce_sweep.

Cases: the 19-knot measured-C-V table (preset), the table scaled by -0.3 / +0.3 dex (the
ENGINE box corners 5e12..5e16 in N_D), the table with its unresolved stretches perturbed
(+0.2 dex in the first 0.5 um, -0.2 dex at 8.5-9.25 um), the preset at 150 / 250 / 400 /
450 K, the table stretched x2 in depth on a 20 um epi (n+ 5 um), and a coarse 3-knot table.

usage: uv run python lab/tabulated_profile/bench/table_parity1d.py
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
from etna.core.device import MEASURED_CV_TABLE, MEASURED_CV_TABULATED_STACK as S
from etna.core.device import donor_table_spec
from etna.core.drift_diffusion import create_dd_device

CV_BIAS = np.array([0.0, -2.0, -5.0, -10.0, -20.0, -40.0, -60.0, -80.0, -100.0])
CCE_BIAS = np.array([0.0, -2.0, -5.0, -10.0, -30.0, -60.0, -100.0])


def scaled(table, dex, lo_um=-1.0, hi_um=1e9):
    return tuple((d, n * 10.0**dex if lo_um <= d <= hi_um else n) for d, n in table)


BASE = dict(
    doping_profile="graded",
    substrate_thickness_cm=S["substrate_thickness_cm"],
    nplus_thickness_cm=S["nplus_thickness_cm"],
    N_D_nplus=S["N_D_nplus"],
    buffer_thickness_cm=None,
    buffer_end_doping=None,
)
T20 = tuple((2 * d, n) for d, n in MEASURED_CV_TABLE)
CASES = [
    ("measured-CV table (preset)", dict(BASE, donor_table=MEASURED_CV_TABLE)),
    (
        "table -0.3 dex (box low)",
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, -0.3)),
    ),
    (
        "table +0.3 dex (box high)",
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, 0.3)),
    ),
    (
        "table +0.2 dex in 0-0.5 um",
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, 0.2, 0.0, 0.5)),
    ),
    (
        "table -0.2 dex in 8.5-9.25 um",
        dict(BASE, donor_table=scaled(MEASURED_CV_TABLE, -0.2, 8.5, 9.25)),
    ),
    ("preset, 150 K", dict(BASE, donor_table=MEASURED_CV_TABLE, T=150.0)),
    ("preset, 250 K", dict(BASE, donor_table=MEASURED_CV_TABLE, T=250.0)),
    ("preset, 400 K", dict(BASE, donor_table=MEASURED_CV_TABLE, T=400.0)),
    ("preset, 450 K", dict(BASE, donor_table=MEASURED_CV_TABLE, T=450.0)),
    (
        "table x2 depth, epi 20 um, n+ 5 um",
        dict(BASE, donor_table=T20, epi_thickness_cm=20e-4, nplus_thickness_cm=5e-4),
    ),
    (
        "3-knot table 1e15 -> 5e13 -> 1e16",
        dict(BASE, donor_table=((0.0, 1e15), (5.0, 5e13), (10.0, 1e16))),
    ),
]


def batch_for(kw):
    epi = kw.get("epi_thickness_cm", 10e-4)
    sub = kw["substrate_thickness_cm"]
    x = node_positions(epi, sub, kw["nplus_thickness_cm"])
    spec = donor_table_spec(epi, kw["donor_table"])
    pb = make_batch(
        x,
        sub,
        kw.get("N_A", 1e19),
        1e15,
        1e14,
        1e-4,
        kw.get("T", 300.0),
        nplus=(sub + epi, kw["N_D_nplus"]),
        buffer=None,
        donor_table=spec,
    )
    return x, pb, epi, sub, spec


def main():
    worst = dict(W=0.0, psi=0.0, cce=0.0, donors=0.0)
    for label, kw in CASES:
        name = f"tpar_{uuid.uuid4().hex[:6]}"
        info = create_dd_device(device_name=name, **kw)
        ref = cv_sweep(info, CV_BIAS)
        g = lambda n: np.array(
            devsim.get_node_model_values(device=name, region="sic", name=n)
        )
        psi_ref, x_ref, don_ref = g("Potential"), g("x"), g("Donors")
        devsim.delete_device(device=name)
        x, pb, epi, sub, spec = batch_for(kw)
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
        )
        nc = len(c["cce_values"])
        dc = np.max(np.abs(ce.cce[:nc, 0] - c["cce_values"]))
        for k, v in (("W", dW), ("psi", dpsi), ("cce", dc), ("donors", ddon)):
            worst[k] = max(worst[k], v)
        print(
            f"RES {label:36s} cv {nref}/{len(CV_BIAS)} trunc {ref['truncated']} conv {res.converged.all()} "
            f"dDonors {ddon:.1e} dW {dW:.2e} dpsi {dpsi:.2e} V | cce {nc}/{len(CCE_BIAS)} conv {ce.converged.all()} "
            f"dCCE {dc:.2e} CCE(-100) {ce.cce[-1, 0]:.6f}",
            flush=True,
        )
    print("RES WORST", {k: f"{v:.2e}" for k, v in worst.items()})


if __name__ == "__main__":
    main()

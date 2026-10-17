"""CCE vs bias of the microdosimeter (tilted 65 deg 5.5 MeV alpha) vs digitised Fig. 4b.

Shape-only comparison: model and data both normalised to their -50 V value.
Variants (sensitivity brackets, NOT tuned to the data): p+ 0.3/1 um, N_D (30 V fit /
1e14), generation shape (Bragg-Kleeman / flat box), metal dead layer (0 / 0.3 um
SiC-eq), SRH lifetimes x0.1/x1/x10, hole mobility from N_A (production) or N_D.

usage: uv run python lab/microdosimeter/bench/cce_bias.py [--out notes/cce_bias.npz]
"""

import itertools
import sys
import time

import numpy as np

sys.path.insert(0, __file__.rsplit("/bench/", 1)[0])
import mdlab  # noqa: E402

N_FIT = 3.0575e14
V = np.array([0, -1, -2, -3, -5, -7.5, -10, -15, -20, -25, -30, -40, -50, -60], float)
PEAK = 1e18  # cm^-3 s^-1, low injection (G*tau_p = 6e11 << N_D)

rows = []
t0 = time.time()
for sub_um in (0.3, 1.0):
    sub = sub_um * 1e-4
    x = mdlab.mesh(sub)
    combos = list(
        itertools.product(
            (N_FIT, 1e14),
            ("bragg", "box"),
            (0.0, 0.3e-4),
            (0.1, 1.0, 10.0),
            ("epi", "p+epi"),
        )
    )
    Gs = []
    for N, shape, td, tau, gen in combos:
        G, lost = mdlab.tilted_generation(
            x, t_dead_cm=td, shape=shape, epi_end=sub + mdlab.EPI_CM
        )
        if (
            gen == "epi"
        ):  # p+ treated as dead layer (run_cce convention; helper default)
            G[x < sub] = 0.0
        Gs.append(G / G.max() * PEAK)
    Gs = np.array(Gs)
    Ns = np.array([c[0] for c in combos])
    pb, _ = mdlab.uniform_batch(sub, Ns)
    for hm in ("N_A", "N_D"):
        cce, conv = mdlab.cce_sweep_custom(
            pb,
            sub,
            V,
            Gs,
            N_D_mob=Ns,
            tau_scale=[c[3] for c in combos],
            hole_mob_from=hm,
        )
        for j, c in enumerate(combos):
            rows.append(
                dict(
                    sub=sub_um,
                    N=c[0],
                    shape=c[1],
                    tdead=c[2] * 1e4,
                    tau=c[3],
                    hmob=hm,
                    gen=c[4],
                    cce=cce[:, j],
                    conv=bool(conv[:, j].all()),
                )
            )
print(f"solved {len(rows)} curves x {len(V)} biases in {time.time() - t0:.1f} s")

iV = [list(V).index(v) for v in mdlab.FIG4B_V]
i50 = list(V).index(-50.0)
data = mdlab.FIG4B_CCE
hdr = "  ".join(f"{v:>6.0f}" for v in mdlab.FIG4B_V)
print(f"{'variant':<44} conv  CCE_abs(-50)  | norm CCE at V = {hdr} | rms dev")
print(
    f"{'DATA (Fig.4b digitised)':<44}       {'':>11}  | "
    + "  ".join(f"{d:6.4f}" for d in data)
)
for r in rows:
    m = r["cce"][iV] / r["cce"][i50]
    rms = np.sqrt(np.mean((m - data) ** 2))
    r["rms"] = rms
    name = f"p+{r['sub']} N{r['N']:.2e} {r['shape']} td{r['tdead']:.1f} tau{r['tau']} mu_p({r['hmob']}) {r['gen']}"
    print(
        f"{name:<50} {str(r['conv'])[0]}     {r['cce'][i50]:.4f}      | "
        + "  ".join(f"{v:6.4f}" for v in m)
        + f" | {rms:.4f}"
    )
np.savez(
    __file__.rsplit("/bench/", 1)[0] + "/notes/cce_bias.npz",
    V=V,
    cce=np.array([r["cce"] for r in rows]),
    labels=np.array(
        [
            f"{r['sub']}|{r['N']:.4e}|{r['shape']}|{r['tdead']}|{r['tau']}|{r['hmob']}|{r['gen']}"
            for r in rows
        ]
    ),
)

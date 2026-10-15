"""Hardware spike: can PyTorch on Apple MPS (Metal) beat vectorized numpy float64
for the inner loop of a BATCHED 1D semiconductor Newton solver?

B independent devices, each a tridiagonal system of N unknowns.

Run with a throwaway venv (see spike_torch_mps.md for exact commands):
    <torch_env>/bin/python spike_torch_mps.py [--quick]

Everything printed is measured on this machine; nothing is estimated.
"""

from __future__ import annotations

import json
import math
import platform
import statistics
import sys
import time

import numpy as np
import torch
import torch.nn.functional as F

QUICK = "--quick" in sys.argv
REPS, WARM = 5, 2
SLOW_SINGLE_S = (
    20.0  # if the first warmup takes longer than this, report that single run
)
B_LIST = [1, 100, 1000, 10000]
N_LIST = [327, 1000]
DENSE_MAX_BYTES = 4e9  # strict: skip dense when matrix alone >= 4 GB (B=1000,N=1000 f32 OOMs on MPS)

OUT_JSON = sys.argv[sys.argv.index("--json") + 1] if "--json" in sys.argv else None

# ----------------------------------------------------------------------------- helpers


def sync_mps():
    torch.mps.synchronize()


def sync_none():
    pass


def bench(fn, sync, reps=REPS, warm=WARM):
    """Median wall time in ms of `reps` runs after `warm` warmups. `sync` is called
    before stopping the clock (torch.mps.synchronize for MPS)."""
    t0 = time.perf_counter()
    fn()
    sync()
    first = time.perf_counter() - t0
    if first > SLOW_SINGLE_S:
        return {
            "median_ms": first * 1e3,
            "min_ms": first * 1e3,
            "max_ms": first * 1e3,
            "runs": 1,
        }
    for _ in range(warm - 1):
        fn()
        sync()
    ts = []
    for _ in range(reps):
        t0 = time.perf_counter()
        fn()
        sync()
        ts.append((time.perf_counter() - t0) * 1e3)
    return {
        "median_ms": statistics.median(ts),
        "min_ms": min(ts),
        "max_ms": max(ts),
        "runs": reps,
    }


def fmt(r):
    if r is None:
        return "skip"
    s = f"{r['median_ms']:.3g}"
    if r.get("runs", REPS) == 1:
        s += "*"
    return s


# ----------------------------------------------------------------------------- problem


def make_mesh(N):
    """N unknowns; N+1 spacings geometrically graded from 3e-8 to 1e-5 (cm)."""
    h = np.geomspace(3e-8, 1e-5, N + 1)
    return h


def make_system(B, N, rng):
    """Return float64 (a, b, c, d) of shape (B, N):
    a = sub-diagonal (a[:,0] unused = 0), c = super-diagonal (c[:,-1] unused = 0),
    b = -(1/h_left + 1/h_right) - d_i, d_i log-uniform in [1e-2, 1e8], rhs random normal.
    """
    h = make_mesh(N)
    inv_l = 1.0 / h[:-1]  # left spacing of unknown i
    inv_r = 1.0 / h[1:]  # right spacing of unknown i
    di = 10.0 ** rng.uniform(-2.0, 8.0, size=(B, N))
    a = np.broadcast_to(inv_l, (B, N)).copy()
    c = np.broadcast_to(inv_r, (B, N)).copy()
    a[:, 0] = 0.0
    c[:, -1] = 0.0
    b = -(inv_l + inv_r)[None, :] - di
    d = rng.standard_normal((B, N))
    return a, b, c, d


def row_scale(a, b, c, d):
    """Divide each row by its diagonal (done in float64 before casting)."""
    return a / b, np.ones_like(b), c / b, d / b


# ----------------------------------------------------------------------------- numpy solvers


def thomas_np(a, b, c, d):
    """Vectorized over B (axis 0), Python loop over N (axis 1). float64."""
    B, N = b.shape
    bp = b.copy()
    dp = d.copy()
    for i in range(1, N):
        w = a[:, i] / bp[:, i - 1]
        bp[:, i] -= w * c[:, i - 1]
        dp[:, i] -= w * dp[:, i - 1]
    x = np.empty_like(d)
    x[:, -1] = dp[:, -1] / bp[:, -1]
    for i in range(N - 2, -1, -1):
        x[:, i] = (dp[:, i] - c[:, i] * x[:, i + 1]) / bp[:, i]
    return x


def pcr_np(a, b, c, d):
    """Parallel cyclic reduction, fully vectorized over (B, N), numpy (extra reference)."""
    N = b.shape[1]
    s = 1
    while s < N:
        pz = ((0, 0), (s, s))
        ap = np.pad(a, pz)
        bp = np.pad(b, pz, constant_values=1.0)
        cp = np.pad(c, pz)
        dp = np.pad(d, pz)
        alpha = -a / bp[:, :N]  # row i-s
        gamma = -c / bp[:, 2 * s :]  # row i+s
        a, c, b, d = (
            alpha * ap[:, :N],
            gamma * cp[:, 2 * s :],
            b + alpha * cp[:, :N] + gamma * ap[:, 2 * s :],
            d + alpha * dp[:, :N] + gamma * dp[:, 2 * s :],
        )
        s *= 2
    return d / b


# ----------------------------------------------------------------------------- torch solvers


def thomas_torch_loop(a, b, c, d):
    """Python loop over N on (B, N) tensors. Uses per-row (B,) views of the (N, B)
    transposed layout so each row is contiguous. Dispatch-bound on MPS."""
    aT, bT, cT, dT = a.T, b.T, c.T, d.T
    N = bT.shape[0]
    bp = [bT[0]]
    dp = [dT[0]]
    for i in range(1, N):
        w = aT[i] / bp[i - 1]
        bp.append(bT[i] - w * cT[i - 1])
        dp.append(dT[i] - w * dp[i - 1])
    x = [None] * N
    x[-1] = dp[-1] / bp[-1]
    for i in range(N - 2, -1, -1):
        x[i] = (dp[i] - cT[i] * x[i + 1]) / bp[i]
    return torch.stack(x, dim=1)


def pcr_torch(a, b, c, d):
    """Parallel cyclic reduction, fully vectorized over (B, N). ceil(log2 N) steps,
    each ~ 4 pads + 12 elementwise kernels."""
    N = b.shape[1]
    s = 1
    while s < N:
        ap = F.pad(a, (s, s))
        bp = F.pad(b, (s, s), value=1.0)
        cp = F.pad(c, (s, s))
        dp = F.pad(d, (s, s))
        alpha = -a / bp[:, :N]
        gamma = -c / bp[:, 2 * s :]
        a, c, b, d = (
            alpha * ap[:, :N],
            gamma * cp[:, 2 * s :],
            b + alpha * cp[:, :N] + gamma * ap[:, 2 * s :],
            d + alpha * dp[:, :N] + gamma * dp[:, 2 * s :],
        )
        s *= 2
    return d / b


def dense_torch(a, b, c, d):
    A = (
        torch.diag_embed(b)
        + torch.diag_embed(a[:, 1:], offset=-1)
        + torch.diag_embed(c[:, :-1], offset=1)
    )
    return torch.linalg.solve(A, d.unsqueeze(-1)).squeeze(-1)


def dense_build(a, b, c):
    return (
        torch.diag_embed(b)
        + torch.diag_embed(a[:, 1:], offset=-1)
        + torch.diag_embed(c[:, :-1], offset=1)
    )


# ----------------------------------------------------------------------------- Newton assembly


def assembly_np(psi, h, dop, s):
    """F and J on interior points (B, N-2). h has N+1 entries (spacings)."""
    # unknown i has left spacing h[i], right spacing h[i+1]; interior i = 1..N-2
    hl = h[1 : N_CUR - 1]
    hr = h[2:N_CUR]
    pm, p0, pp = psi[:, :-2], psi[:, 1:-1], psi[:, 2:]
    em, ep = np.exp(-p0), np.exp(p0)
    Fv = (pp - p0) / hr - (p0 - pm) / hl + s * (em - ep + dop[:, 1:-1])
    Jv = -(1.0 / hr + 1.0 / hl) - s * (em + ep)
    return Fv, Jv


def assembly_torch(psi, hl, hr, dop_int, s):
    pm, p0, pp = psi[:, :-2], psi[:, 1:-1], psi[:, 2:]
    em, ep = torch.exp(-p0), torch.exp(p0)
    Fv = (pp - p0) / hr - (p0 - pm) / hl + s * (em - ep + dop_int)
    Jv = -(1.0 / hr + 1.0 / hl) - s * (em + ep)
    return Fv, Jv


N_CUR = 327  # set per run for assembly_np slicing


# ----------------------------------------------------------------------------- full Newton step


def newton_step_np(psi, h, dop, s):
    Fv, Jv = assembly_np(psi, h, dop, s)
    B, N = psi.shape
    a = np.zeros((B, N))
    c = np.zeros((B, N))
    b = np.ones((B, N))
    rhs = np.zeros((B, N))
    a[:, 1:-1] = 1.0 / h[1 : N - 1]
    c[:, 1:-1] = 1.0 / h[2:N]
    b[:, 1:-1] = Jv
    rhs[:, 1:-1] = -Fv
    dx = thomas_np(a, b, c, rhs)
    return psi + dx


def newton_step_torch(psi, hl, hr, dop_int, s, solver):
    Fv, Jv = assembly_torch(psi, hl, hr, dop_int, s)
    B, N = psi.shape
    a = F.pad(hl.reciprocal().expand(B, N - 2), (1, 1))
    c = F.pad(hr.reciprocal().expand(B, N - 2), (1, 1))
    b = F.pad(Jv, (1, 1), value=1.0)
    rhs = F.pad(-Fv, (1, 1))
    dx = solver(a, b, c, rhs)
    return psi + dx


# ----------------------------------------------------------------------------- main

results = {
    "env": {},
    "float64_mps": None,
    "tridiag": {},
    "accuracy": {},
    "assembly": {},
    "newton": {},
}


def main():
    global N_CUR
    env = results["env"]
    env["python"] = sys.version.split()[0]
    env["torch"] = torch.__version__
    env["numpy"] = np.__version__
    env["platform"] = platform.platform()
    env["machine"] = platform.machine()
    env["mps_available"] = torch.backends.mps.is_available()
    env["mps_built"] = torch.backends.mps.is_built()
    env["threads"] = torch.get_num_threads()
    print("ENV", json.dumps(env, indent=1))
    if not env["mps_available"]:
        print("MPS not available; abort")
        return

    mps = torch.device("mps")

    # 1. float64 on MPS
    try:
        t = torch.zeros(4, dtype=torch.float64, device=mps)
        results["float64_mps"] = {"supported": True, "msg": str(t.dtype)}
    except Exception as e:  # noqa: BLE001
        results["float64_mps"] = {
            "supported": False,
            "error": f"{type(e).__name__}: {e}",
        }
    print("FLOAT64_MPS", json.dumps(results["float64_mps"], indent=1))

    rng = np.random.default_rng(0)
    b_list = [1, 100, 1000] if QUICK else B_LIST

    # 2-5. tridiagonal solvers
    for N in N_LIST:
        h = make_mesh(N)
        print(
            f"\n=== N={N}: h in [{h.min():.3g}, {h.max():.3g}] cm, total length {h.sum():.4g} cm ==="
        )
        for B in b_list:
            key = f"B={B},N={N}"
            a, b, c, d = make_system(B, N, rng)
            x64 = thomas_np(a, b, c, d)  # reference
            r = {}

            # numpy float64 Thomas
            r["np_f64_thomas"] = bench(lambda: thomas_np(a, b, c, d), sync_none)
            # numpy float64 PCR (extra)
            r["np_f64_pcr"] = bench(lambda: pcr_np(a, b, c, d), sync_none)

            # torch MPS float32 tensors
            ta, tb, tc, td = (
                torch.tensor(v, dtype=torch.float32, device=mps) for v in (a, b, c, d)
            )
            sync_mps()

            # Thomas loop on MPS (dispatch bound): cap cost at large N*B -- still run all
            r["mps_f32_thomas_loop"] = bench(
                lambda: thomas_torch_loop(ta, tb, tc, td), sync_mps
            )
            r["mps_f32_pcr"] = bench(lambda: pcr_torch(ta, tb, tc, td), sync_mps)

            # PCR compiled
            try:
                pcr_c = torch.compile(pcr_torch)
                t0 = time.perf_counter()
                pcr_c(ta, tb, tc, td)
                sync_mps()
                compile_s = time.perf_counter() - t0
                r["mps_f32_pcr_compiled"] = bench(
                    lambda: pcr_c(ta, tb, tc, td), sync_mps
                )
                r["mps_f32_pcr_compiled"]["compile_s"] = compile_s
            except Exception as e:  # noqa: BLE001
                r["mps_f32_pcr_compiled"] = None
                r["mps_f32_pcr_compiled_error"] = f"{type(e).__name__}: {str(e)[:300]}"

            # dense MPS float32 (build + solve)
            dense_bytes = B * N * N * 4
            if dense_bytes < DENSE_MAX_BYTES and not (B == 10000):
                r["mps_f32_dense"] = bench(
                    lambda: dense_torch(ta, tb, tc, td), sync_mps
                )
                A = dense_build(ta, tb, tc)
                sync_mps()
                r["mps_f32_dense_solve_only"] = bench(
                    lambda: torch.linalg.solve(A, td.unsqueeze(-1)), sync_mps
                )
                del A
            else:
                r["mps_f32_dense"] = None
                r["mps_f32_dense_solve_only"] = None
            # dense CPU float64 (build + solve)
            if B * N * N * 8 < DENSE_MAX_BYTES and not (B == 10000):
                ca, cb, cc, cd = (
                    torch.tensor(v, dtype=torch.float64) for v in (a, b, c, d)
                )
                r["cpu_f64_dense"] = bench(
                    lambda: dense_torch(ca, cb, cc, cd), sync_none
                )
                del ca, cb, cc, cd
            else:
                r["cpu_f64_dense"] = None

            results["tridiag"][key] = r
            print(
                f"{key:>14}  np64_thomas {fmt(r['np_f64_thomas']):>8}  np64_pcr {fmt(r['np_f64_pcr']):>8}  "
                f"mps_thomas {fmt(r['mps_f32_thomas_loop']):>8}  mps_pcr {fmt(r['mps_f32_pcr']):>8}  "
                f"mps_pcr_c {fmt(r['mps_f32_pcr_compiled']):>8}  mps_dense {fmt(r['mps_f32_dense']):>8}  "
                f"(solve-only {fmt(r['mps_f32_dense_solve_only'])})  cpu64_dense {fmt(r['cpu_f64_dense']):>8}",
                flush=True,
            )

            # accuracy vs numpy float64 (max over systems of relative inf-norm error)
            def relerr(x32):
                x32 = np.asarray(x32, dtype=np.float64)
                num = np.abs(x32 - x64).max(axis=1)
                den = np.abs(x64).max(axis=1)
                return float((num / den).max())

            acc = {}
            for scaled in (False, True):
                aa, bb, cc2, dd = row_scale(a, b, c, d) if scaled else (a, b, c, d)
                sa, sb, sc, sd = (
                    torch.tensor(v, dtype=torch.float32, device=mps)
                    for v in (aa, bb, cc2, dd)
                )
                tag = "scaled" if scaled else "raw"
                acc[f"np_f32_thomas_{tag}"] = relerr(
                    thomas_np(*(v.astype(np.float32) for v in (aa, bb, cc2, dd)))
                )
                acc[f"mps_f32_thomas_{tag}"] = relerr(
                    thomas_torch_loop(sa, sb, sc, sd).cpu().numpy()
                )
                acc[f"mps_f32_pcr_{tag}"] = relerr(
                    pcr_torch(sa, sb, sc, sd).cpu().numpy()
                )
                acc[f"np_f64_pcr_{tag}"] = relerr(pcr_np(aa, bb, cc2, dd))
                if r["mps_f32_dense"] is not None:
                    acc[f"mps_f32_dense_{tag}"] = relerr(
                        dense_torch(sa, sb, sc, sd).cpu().numpy()
                    )
                del sa, sb, sc, sd
            results["accuracy"][key] = acc
            print("   accuracy:", {k: f"{v:.2e}" for k, v in acc.items()}, flush=True)
            del ta, tb, tc, td
            torch.mps.empty_cache()

    # 6. elementwise assembly
    print("\n=== Newton assembly (elementwise on (B, N)) ===")
    s_val = 1e-3
    for N in N_LIST:
        N_CUR = N
        h = make_mesh(N)
        hl_np = h[1 : N - 1]
        hr_np = h[2:N]
        for B in b_list:
            key = f"B={B},N={N}"
            psi = rng.uniform(-20, 20, size=(B, N))
            dop = rng.uniform(-1e3, 1e3, size=(B, N))
            r = {}
            r["np_f64"] = bench(lambda: assembly_np(psi, h, dop, s_val), sync_none)
            tpsi = torch.tensor(psi, dtype=torch.float32, device=mps)
            tdop = torch.tensor(dop[:, 1:-1], dtype=torch.float32, device=mps)
            thl = torch.tensor(hl_np, dtype=torch.float32, device=mps)
            thr = torch.tensor(hr_np, dtype=torch.float32, device=mps)
            sync_mps()
            r["mps_f32_eager"] = bench(
                lambda: assembly_torch(tpsi, thl, thr, tdop, s_val), sync_mps
            )
            try:
                asm_c = torch.compile(assembly_torch)
                t0 = time.perf_counter()
                asm_c(tpsi, thl, thr, tdop, s_val)
                sync_mps()
                compile_s = time.perf_counter() - t0
                r["mps_f32_compiled"] = bench(
                    lambda: asm_c(tpsi, thl, thr, tdop, s_val), sync_mps
                )
                r["mps_f32_compiled"]["compile_s"] = compile_s
            except Exception as e:  # noqa: BLE001
                r["mps_f32_compiled"] = None
                r["mps_f32_compiled_error"] = f"{type(e).__name__}: {str(e)[:300]}"
            # torch CPU float64 eager (extra reference)
            cpsi = torch.tensor(psi, dtype=torch.float64)
            cdop = torch.tensor(dop[:, 1:-1], dtype=torch.float64)
            chl = torch.tensor(hl_np, dtype=torch.float64)
            chr_ = torch.tensor(hr_np, dtype=torch.float64)
            r["cpu_f64_torch_eager"] = bench(
                lambda: assembly_torch(cpsi, chl, chr_, cdop, s_val), sync_none
            )
            results["assembly"][key] = r
            print(
                f"{key:>14}  np64 {fmt(r['np_f64']):>8}  mps32_eager {fmt(r['mps_f32_eager']):>8}  "
                f"mps32_compiled {fmt(r['mps_f32_compiled']):>8}  torch_cpu64 {fmt(r['cpu_f64_torch_eager']):>8}",
                flush=True,
            )
            del tpsi, tdop
            torch.mps.empty_cache()

    # 7. full Newton iteration at N=327, B in {1000, 10000}
    print("\n=== Full Newton iteration (assembly + tridiag solve + update), N=327 ===")
    N = 327
    N_CUR = N
    h = make_mesh(N)
    hl_np = h[1 : N - 1]
    hr_np = h[2:N]
    for B in [1000] if QUICK else [1000, 10000]:
        key = f"B={B},N={N}"
        psi = rng.uniform(-20, 20, size=(B, N))
        dop = rng.uniform(-1e3, 1e3, size=(B, N))
        r = {}
        r["np_f64"] = bench(lambda: newton_step_np(psi, h, dop, s_val), sync_none)
        tpsi = torch.tensor(psi, dtype=torch.float32, device=mps)
        tdop = torch.tensor(dop[:, 1:-1], dtype=torch.float32, device=mps)
        thl = torch.tensor(hl_np, dtype=torch.float32, device=mps)
        thr = torch.tensor(hr_np, dtype=torch.float32, device=mps)
        sync_mps()
        r["mps_f32_pcr_eager"] = bench(
            lambda: newton_step_torch(tpsi, thl, thr, tdop, s_val, pcr_torch), sync_mps
        )
        try:
            step_c = torch.compile(
                lambda p: newton_step_torch(p, thl, thr, tdop, s_val, pcr_torch)
            )
            t0 = time.perf_counter()
            step_c(tpsi)
            sync_mps()
            compile_s = time.perf_counter() - t0
            r["mps_f32_pcr_compiled"] = bench(lambda: step_c(tpsi), sync_mps)
            r["mps_f32_pcr_compiled"]["compile_s"] = compile_s
        except Exception as e:  # noqa: BLE001
            r["mps_f32_pcr_compiled"] = None
            r["mps_f32_pcr_compiled_error"] = f"{type(e).__name__}: {str(e)[:300]}"
        if B * N * N * 4 < DENSE_MAX_BYTES and B != 10000:
            r["mps_f32_dense"] = bench(
                lambda: newton_step_torch(tpsi, thl, thr, tdop, s_val, dense_torch),
                sync_mps,
            )
        else:
            r["mps_f32_dense"] = None
        # accuracy of the torch step vs numpy step (one iteration from the same psi)
        psi_np = newton_step_np(psi, h, dop, s_val)
        psi_t = (
            newton_step_torch(tpsi, thl, thr, tdop, s_val, pcr_torch)
            .cpu()
            .numpy()
            .astype(np.float64)
        )
        dx_np = psi_np - psi
        dx_t = psi_t - psi
        r["dx_relerr_pcr"] = float(
            (np.abs(dx_t - dx_np).max(axis=1) / np.abs(dx_np).max(axis=1)).max()
        )
        results["newton"][key] = r
        print(
            f"{key:>14}  np64 {fmt(r['np_f64']):>8}  mps32_pcr_eager {fmt(r['mps_f32_pcr_eager']):>8}  "
            f"mps32_pcr_compiled {fmt(r['mps_f32_pcr_compiled']):>8}  mps32_dense {fmt(r['mps_f32_dense']):>8}  "
            f"dx_relerr(pcr) {r['dx_relerr_pcr']:.2e}",
            flush=True,
        )
        del tpsi, tdop
        torch.mps.empty_cache()

    if OUT_JSON:
        with open(OUT_JSON, "w") as f:
            json.dump(results, f, indent=1, default=str)
        print("wrote", OUT_JSON)


if __name__ == "__main__":
    main()

"""Batched 1D engine: the fused carrier-solve kernel equals the numpy assembly + Thomas
reference, and the serial / prange drivers (chosen by batch size) give identical results."""

import numpy as np
import pytest

import etna.core.batched.poisson as po
import etna.core.batched.transport as tr

N = 300


def _batch(B):
    x = np.linspace(0.0, 11e-4, N)
    don = np.where(x > 0.3e-4, 3.4e13, 0.0)[None, :].repeat(B, 0)
    acc = np.where(x <= 0.3e-4, 1e19, 0.0)[None, :].repeat(B, 0)
    pb = po.PoissonBatch(x, don, acc, 5e-9, 0.02585)
    V = np.linspace(2.0, 20.0, B)
    res = pb.solve(pb.neutral_guess(0.0, V), V_cathode=V)
    return pb, res


@pytest.mark.parametrize("B", [1, 48])
@pytest.mark.parametrize("auger", [False, True])
def test_fused_carrier_solve_equals_numpy_reference(B, auger):
    pb, res = _batch(B)
    C = (2.8e-31, 9.9e-32) if auger else (0.0, 0.0)
    cb = tr.ContinuityBatch(pb, 1000.0, 100.0, 1e-7, 1e-7, C_n=C[0], C_p=C[1])
    psi_T = res.psi.T.copy()
    n_T, p_T = res.n.T.copy(), res.p.T.copy()
    G = np.full((N, B), 1e18)
    denom = 1e-7 * (n_T + pb.n_i) + 1e-7 * (p_T + pb.n_i)
    aug = cb.C_n[None, :] * n_T + cb.C_p[None, :] * p_T if auger else None
    for which, other, bc in (("n", p_T, cb.n_bc), ("p", n_T, cb.p_bc)):
        x = cb._solve_carrier(psi_T, other, which, G, denom, bc, aug)
        r = cb._solve_carrier_numpy(psi_T, other, which, G, denom, bc, aug)
        np.testing.assert_array_equal(x, r)


def test_serial_and_prange_drivers_identical():
    B = 48
    assert po.use_prange(B, N, po.NEWTON_PAR_NODES)
    assert not po.use_prange(1, N, po.NEWTON_PAR_NODES)
    pb, res = _batch(B)
    u0 = np.ascontiguousarray(res.psi / pb.V_t[:, None] * (1 + 1e-6))
    un = np.full((B, N), 5.0 / 0.02585)
    up = np.zeros((B, N))
    out = []
    for f in (po._newton_seq, po._newton_par):
        u = u0.copy()
        it = np.zeros(B, np.int64)
        cv = np.zeros(B, np.bool_)
        f(
            u,
            un,
            up,
            pb._D_hat,
            pb.log_ni_hat,
            pb.lam2,
            pb._inv_h,
            pb.v,
            1e-10,
            50,
            20.0,
            it,
            cv,
        )
        out.append((u, it, cv))
    for a, b in zip(*out):
        np.testing.assert_array_equal(a, b)
    nh = np.ascontiguousarray(res.n / pb.N0[:, None])
    ph = np.ascontiguousarray(res.p / pb.N0[:, None])
    c1, c2 = np.empty(B), np.empty(B)
    po._dpsi_dvc_seq(nh, ph, pb.lam2, pb._inv_h, pb.v, c1)
    po._dpsi_dvc_par(nh, ph, pb.lam2, pb._inv_h, pb.v, c2)
    np.testing.assert_array_equal(c1, c2)

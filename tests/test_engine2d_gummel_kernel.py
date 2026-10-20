"""Compiled illuminated Gummel loop of the 2D engine: column-oriented band solve and the
parallel stripe batch reproduce the row-oriented reference (the former Python loop's
expression order) to rounding."""

import math

import numba
import numpy as np
import pytest

from etna.api.device import DeviceConfig, device_kwargs

from etna.core.engine2d.banded import band_lu, band_solve, band_solve_cols, band_split
from etna.core.engine2d.device import Device2D
from etna.core.engine2d.grid import grid_from_devsim
from etna.core.engine2d.kernels import cathode_current_balanced
from etna.core.engine2d.scan import _generation, _uniform_epi_generation


def test_band_solve_cols_matches_band_solve():
    rng = np.random.default_rng(1)
    N, w = 300, 7
    A = rng.uniform(-1.0, 0.0, (N, 2 * w + 1))
    A[:, w] = 2 * w + 1.0 + rng.uniform(0, 1, N)  # diagonally dominant
    for k in range(N):  # outside the matrix
        for m in range(2 * w + 1):
            if not 0 <= k + m - w < N:
                A[k, m] = 0.0
    band_lu(A, w)
    b = rng.standard_normal(N)
    x_ref = band_solve(A, w, b)
    x = b.copy()
    band_solve_cols(*band_split(A, w), w, x)
    assert np.max(np.abs(x - x_ref) / np.abs(x_ref)) < 1e-13


def _reference_illuminated(dev, G, tol=1e-10, tol_I=1e-9, max_iter=60):
    """The former Python Gummel loop of Device2D.illuminated_current (row band_solve)."""
    g = dev.g2
    w, lnni = g.nx, dev.lnni
    u, n, p = dev.u.copy(), dev.n.copy(), dev.p.copy()
    R, F, dummy = np.empty(g.N), np.empty(g.N), np.empty((1, 1))
    I_prev = None
    for _ in range(max_iter):
        dev._cont(0, u, n, p, G, dummy, R, False)
        n += band_solve(dev.LU_n, w, -R)
        dev._cont(1, u, n, p, G, dummy, R, False)
        p += band_solve(dev.LU_p, w, -R)
        np.maximum(n, 1e-300, out=n)
        np.maximum(p, 1e-300, out=p)
        un = u - (np.log(n) - lnni)
        up = u + (np.log(p) - lnni)
        dev._poisson(
            u, un, up, lnni, dev.lam, g.hx, g.hy, g.dxd, g.dyd, g.nx, g.ny, dummy, F,
            False,
        )  # fmt: skip
        du = band_solve(dev.LU_psi, w, -F)
        u += du
        e = np.exp(du)
        n *= e
        p /= e
        Ie, Ih = cathode_current_balanced(
            u, n, p, G, dev.ni, dev.taun, dev.taup, dev.kn_v, dev.kp_v, g.hy, g.dxd,
            g.dyd, g.nx, g.ny,
        )  # fmt: skip
        I = Ie + Ih
        dI = abs(I - I_prev) / abs(I) if I_prev is not None and I != 0 else math.inf
        I_prev = I
        if np.abs(du).max() < tol and dI < tol_I:
            return I
    raise AssertionError("reference Gummel did not converge")


@pytest.fixture(scope="module", params=[False, True], ids=["static", "ionization"])
def dark_device(request):
    kw = {}
    if request.param:
        kw["incomplete_ionization"] = device_kwargs(
            DeviceConfig(incomplete_ionization=True)
        )["incomplete_ionization"]
    grid, info = grid_from_devsim(half_width_um=25.0, epi_thickness_cm=10e-4, **kw)
    dev = Device2D(grid, info)
    assert (dev.ion is not None) == request.param
    dev.dark(20.0)
    dev.factor_continuity()
    xs = np.linspace(0, info["half_width_cm"], 3)
    Gs = [_generation(dev, x, 1e18, 2e-4) for x in xs]
    Gs.append(_uniform_epi_generation(dev, 1e18))
    return dev, np.array(Gs)


def test_column_solve_and_batch_match_reference(dark_device, monkeypatch):
    dev, Gs = dark_device
    legacy = np.array([_reference_illuminated(dev, G) for G in Gs])
    dev.column_solve = False
    dev.gummel_iterations.clear()
    ref = np.array([dev.illuminated_current(G) for G in Gs])
    ref_its = list(dev.gummel_iterations)
    # the compiled loop with the row solve is the former Python loop, bit for bit
    np.testing.assert_array_equal(ref, legacy)
    dev.column_solve = True
    dev.gummel_iterations.clear()
    col = np.array([dev.illuminated_current(G) for G in Gs])
    dev.gummel_iterations.clear()
    # xdist workers pin NUMBA_NUM_THREADS=1: force the thread-pool branch
    monkeypatch.setattr(numba, "get_num_threads", lambda: 4)
    bat = dev.illuminated_currents(Gs)
    assert dev.gummel_iterations == ref_its
    assert np.max(np.abs(col / ref - 1)) <= 1e-13
    # the batch runs the same per-generation arithmetic as single calls
    np.testing.assert_array_equal(bat, col)
    assert np.isfinite(dev.last_I_direct)


def test_not_converged_raises(dark_device):
    from etna.core.engine2d.device import EngineNotConverged

    dev, Gs = dark_device
    with pytest.raises(EngineNotConverged):
        dev.illuminated_current(Gs[0], max_iter=1)
    with pytest.raises(EngineNotConverged):
        dev.illuminated_currents(Gs, max_iter=1)

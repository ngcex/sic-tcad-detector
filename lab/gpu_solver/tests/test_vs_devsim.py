"""Validation of the batched solvers against ETNA/devsim on identical meshes.

Tolerances are explicit; a failure here means the solver is NOT done.
Requires the ``ref`` extra (devsim). Run: uv run pytest -q
"""

import numpy as np
import pytest

from gpusolver.extract import depletion_width, electric_field
from gpusolver.physics import graded_doping, ionized_acceptor_fraction
from gpusolver.poisson_np import PoissonBatch
from gpusolver.ref_devsim import build_ref_device, delete_ref_device, ref_cv_sweep
from gpusolver import physics
from gpusolver.sweep import cce_sweep, make_batch, reverse_sweep

V_PTS = np.arange(0.0, -61.0, -5.0)

# devsim converges on relative update 1e-10; float64 agreement is at rounding level.
TOL_PSI_F64 = 1e-10  # V, max abs
TOL_W_F64 = 1e-10  # relative
TOL_E_F64 = 1e-10  # relative to E_max
# float32 GPU path: potential resolution ~ |psi|*eps32 ~ 60 V * 1.2e-7 ~ 1e-5 V
TOL_PSI_F32 = 2e-4  # V
TOL_W_F32 = 1e-4  # relative
TOL_E_F32 = 1e-4  # relative to E_max


@pytest.fixture(scope="module")
def ref():
    r = build_ref_device()
    r.cv = ref_cv_sweep(r, V_PTS)  # DD solver, 0.5 V internal ramp, ends at -60 V
    yield r
    delete_ref_device(r)


def test_doping_profile_matches_devsim(ref):
    d, a = graded_doping(
        ref.x, ref.junction_pos, ref.acceptors.max(), 2.93e15, 8.82e13, 0.987e-4
    )
    assert np.allclose(d[0], ref.donors, rtol=1e-12, atol=0)
    assert np.array_equal(a[0], ref.acceptors)
    assert np.isclose(
        1e19 * ionized_acceptor_fraction(1e19, 300.0), ref.acceptors.max(), rtol=1e-12
    )


def test_equilibrium_potential_numpy(ref):
    pb = PoissonBatch(ref.x, ref.donors[None], ref.acceptors[None], ref.n_i, ref.V_t)
    res = pb.solve(pb.neutral_guess(), tol=1e-8)
    assert res.converged.all()
    assert np.abs(res.psi[0] - ref.psi_eq).max() < TOL_PSI_F64


def test_reverse_sweep_numpy_matches_dd(ref):
    """W(V), psi(x) and E(x) at -60 V vs devsim drift-diffusion cv_sweep."""
    assert not ref.cv["truncated"]
    pb = make_batch(ref.x, ref.junction_pos, 1e19, 2.93e15, 8.82e13, 0.987e-4)
    sw = reverse_sweep(pb, V_PTS, ref.junction_pos, ref.epi_thickness_cm)
    assert sw.converged.all()
    W_ref = ref.cv["depletion_widths"]
    assert np.abs(sw.W[:, 0] - W_ref).max() / W_ref.min() < TOL_W_F64
    assert np.abs(sw.psi[-1, 0] - ref.cv["Potential"]).max() < TOL_PSI_F64
    E_ref = -np.diff(ref.cv["Potential"]) / np.diff(ref.x)
    E = electric_field(ref.x, sw.psi[-1])[0]
    assert np.abs(E - E_ref).max() / np.abs(E_ref).max() < TOL_E_F64


def test_reverse_sweep_numba_matches_dd(ref):
    pb = make_batch(ref.x, ref.junction_pos, 1e19, 2.93e15, 8.82e13, 0.987e-4, backend="numba")
    sw = reverse_sweep(pb, V_PTS, ref.junction_pos, ref.epi_thickness_cm)
    assert sw.converged.all()
    W_ref = ref.cv["depletion_widths"]
    assert np.abs(sw.W[:, 0] - W_ref).max() / W_ref.min() < TOL_W_F64
    assert np.abs(sw.psi[-1, 0] - ref.cv["Potential"]).max() < TOL_PSI_F64


def test_cce_numba_matches_etna(ref_cce):
    r = _cce("numba", 1e-8)
    assert r.converged.all()
    assert np.abs(r.cce[:, 0] - ref_cce["cce_values"]).max() < TOL_CCE_F64


def test_reverse_sweep_batch_broadcast(ref):
    """A batch of 3 dopings: the middle one must equal the single-device result."""
    nb = np.array([5e13, 8.82e13, 2e14])
    pb = make_batch(ref.x, ref.junction_pos, 1e19, 2.93e15, nb, 0.987e-4)
    sw = reverse_sweep(pb, V_PTS, ref.junction_pos, ref.epi_thickness_cm)
    assert sw.converged.all()
    W_ref = ref.cv["depletion_widths"]
    assert np.abs(sw.W[:, 1] - W_ref).max() / W_ref.min() < TOL_W_F64
    # physics sanity: lower bulk doping -> wider depletion at moderate bias
    assert sw.W[1, 0] > sw.W[1, 1] > sw.W[1, 2]


@pytest.mark.mlx
def test_reverse_sweep_mlx_float32(ref):
    pytest.importorskip("mlx.core")
    pb = make_batch(
        ref.x, ref.junction_pos, 1e19, 2.93e15, 8.82e13, 0.987e-4, backend="mlx"
    )
    sw = reverse_sweep(pb, V_PTS, ref.junction_pos, ref.epi_thickness_cm, tol=1e-4)
    assert sw.converged.all()
    W_ref = ref.cv["depletion_widths"]
    assert np.abs(sw.W[:, 0] - W_ref).max() / W_ref.min() < TOL_W_F32
    assert np.abs(sw.psi[-1, 0] - ref.cv["Potential"]).max() < TOL_PSI_F32
    E_ref = -np.diff(ref.cv["Potential"]) / np.diff(ref.x)
    E = electric_field(ref.x, sw.psi[-1])[0]
    assert np.abs(E - E_ref).max() / np.abs(E_ref).max() < TOL_E_F32


CCE_V = np.array([0.0, -5.0, -10.0, -40.0])
TOL_CCE_F64 = 1e-6  # absolute, CCE in [0, 1]
TOL_CCE_F32 = 1e-5


@pytest.fixture(scope="module")
def ref_cce():
    from etna.core.charge_collection import cce_vs_bias

    from gpusolver.ref_devsim import quiet

    with quiet():
        r = cce_vs_bias(CCE_V, device_kwargs=dict(N_D_junction=2.93e15, N_D_bulk=8.82e13, L_transition=0.987e-4))
    assert not r["truncated"]
    return r


def _cce(backend, tol):
    X = np.load(__file__.replace("tests/test_vs_devsim.py", "bench/mesh_devsim_default.npy"))
    pb = make_batch(X, 1e-4, 1e19, 2.93e15, 8.82e13, 0.987e-4, backend=backend)
    G = physics.alpha_generation(X, 1e-4)
    mu_n, mu_p = physics.mobility(1.07e15, 300.0, "electron"), physics.mobility(1e19, 300.0, "hole")
    return cce_sweep(pb, CCE_V, 1e-4, G, mu_n, mu_p, 1e-9, 6e-7, tol=tol)


def test_cce_numpy_matches_etna(ref_cce):
    r = _cce("numpy", 1e-8)
    assert r.converged.all()
    assert np.isclose(r.J_generated, ref_cce["I_generated"], rtol=1e-12)
    assert np.abs(r.cce[:, 0] - ref_cce["cce_values"]).max() < TOL_CCE_F64


@pytest.mark.mlx
def test_cce_mlx_matches_etna(ref_cce):
    pytest.importorskip("mlx.core")
    r = _cce("mlx", 1e-4)
    assert r.converged.all()
    assert np.abs(r.cce[:, 0] - ref_cce["cce_values"]).max() < TOL_CCE_F32

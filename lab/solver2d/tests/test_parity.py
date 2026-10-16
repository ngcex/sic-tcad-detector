"""Parity of the 2D engine vs devsim dumps (notes/*.npz, refs ref_*_50.npz in notes/).

Run from the repo root: uv run python -m pytest lab/solver2d/tests -q
(needs notes/mesh_hw{25,50}.npz from bench/dump_mesh.py and notes/scan_hw50_V50.npz
from bench/devsim_scan_dump.py 50 50 10).
"""

import os
import sys
from pathlib import Path

import numpy as np
import pytest

HERE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(HERE))
from solver2d.engine import Device2D  # noqa: E402
from solver2d.grid import grid_from_nodes  # noqa: E402

NOTES = HERE / "notes"
REF = Path(os.environ.get("SOLVER2D_REF", NOTES))  # ref_{50,150}_50.npz are in notes/


def _grid(hw):
    d = np.load(NOTES / f"mesh_hw{hw}.npz")
    return d, grid_from_nodes(d["x"], d["y"], d["elements"])


@pytest.mark.parametrize("hw", [25, 50])
def test_equilibrium(hw):
    d, g = _grid(hw)
    dev = Device2D(g)
    D_k = np.repeat(dev.D, g.nx)
    assert (
        np.max(
            np.abs(g.to_devsim(D_k) - d["NetDoping"]) / np.abs(d["NetDoping"]).clip(1)
        )
        < 1e-12
    )
    psi2d, _ = dev.equilibrium_2d()
    assert np.abs(g.to_devsim(psi2d) - d["Potential"]).max() < 1e-6
    psi = dev.dark(0.0)
    assert np.abs(g.to_devsim(psi) - d["Potential"]).max() < 1e-6


def _carrier_err(a, b):
    m = (b > 1.0) | (a > 1.0)
    return np.abs(np.log10(np.maximum(a[m], 1e-300)) - np.log10(b[m])).max()


@pytest.mark.parametrize("hw", [50, 150])
def test_dark_50V(hw):
    ref = np.load(REF / f"ref_{hw}_50.npz")
    from solver2d.grid import grid_from_devsim

    g = grid_from_devsim(float(hw))
    dev = Device2D(g)
    psi = dev.dark(50.0)
    assert np.abs(g.to_devsim(psi) - ref["Potential"]).max() < 1e-6
    assert _carrier_err(g.to_devsim(dev.n), ref["Electrons"]) < 1e-6
    assert _carrier_err(g.to_devsim(dev.p), ref["Holes"]) < 1e-6


def test_cce_scan_hw50():
    """10-point lateral scan vs devsim (notes/scan_hw50_V50.npz, bench/devsim_scan_dump.py 50 50 10)."""
    d = np.load(NOTES / "scan_hw50_V50.npz")
    _, g = _grid(50)
    dev = Device2D(g)
    dev.dark(50.0)
    dev.factor_continuity()
    _, cce = dev.cce_lateral_scan(10)
    assert np.abs(cce - d["cce"]).max() < 1e-5

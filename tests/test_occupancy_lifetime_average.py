"""Issue #13 (audit v8 M2): the default-path occupancy lifetimes are a LENGTH mean over
the epi (junction..epi_end), not a node-count mean over epi + n+; the CCE no longer
depends on the mesh grading (.planning/AUDIT_v8_MODELS_DAMAGE.md)."""

from __future__ import annotations

import pytest

from etna.api.device import DeviceConfig


# ------------------------------------------------------------------ issue #13 (M2)
def test_occupancy_lifetimes_are_a_length_mean():
    """Same linear availability profile on a graded and a uniform mesh: the length mean
    is mesh independent (trapezoid exact for a linear profile), the node mean is not."""
    import numpy as np

    from etna.core.radiation_traps import (
        DEFAULT_TRAP_MODEL,
        occupancy_lifetimes,
        trap_fields,
    )

    L = 50e-4
    graded = L * np.linspace(0.0, 1.0, 201) ** 3  # dense near x = 0
    uniform = np.linspace(0.0, L, 51)
    out = {}
    for name, x in (("graded", graded), ("uniform", uniform)):
        mask = np.ones(x.size, dtype=bool)
        fl = trap_fields(DEFAULT_TRAP_MODEL, 1e13, mask, 300.0, 5e-9)
        f = np.broadcast_to(x / L, fl.N_t.shape).copy()  # electron occupancy 0 -> 1
        out[name] = occupancy_lifetimes(
            DEFAULT_TRAP_MODEL, fl, f, mask, 300.0, 1e-6, 1e-6, x=x
        )
        out[name + "_node"] = occupancy_lifetimes(
            DEFAULT_TRAP_MODEL, fl, f, mask, 300.0, 1e-6, 1e-6
        )
    assert out["graded"] == pytest.approx(out["uniform"], rel=1e-12)
    # the plain node mean (the pre-fix default) moves with the grading
    assert abs(out["graded_node"][0] / out["uniform_node"][0] - 1.0) > 0.2


def test_default_cce_is_independent_of_the_junction_mesh(monkeypatch):
    """Default path ("auto" -> occupancy lifetimes), 62 MeV, -40 V, 1e14 p/cm^2: four
    times finer spacing at the 1-nm mesh lines (junction, n-/n+ step) moved the CCE by
    1.0e-4 with the node mean (0.9922168 -> 0.9921152), by 6e-7 with the epi length
    mean; the legacy-lifetime control moves by 7e-7 (lab/audit_v8/m2_mesh_refinement.py).
    """
    import devsim
    import numpy as np

    from etna import run_radiation_damage

    def run():
        r = run_radiation_damage(
            DeviceConfig(), fluences=np.array([1e14]), V_bias=-40.0,
            proton_energy_MeV=62.0,
        )  # fmt: skip
        assert r.metadata["damage_model"] == "charged_traps"
        return float(r.y[0])

    base = run()
    orig = devsim.add_1d_mesh_line

    def finer(**kw):
        if kw["ps"] <= 1e-6:
            kw = dict(kw, ps=kw["ps"] * 0.25)
        return orig(**kw)

    monkeypatch.setattr(devsim, "add_1d_mesh_line", finer)
    fine = run()
    assert abs(fine - base) < 5e-6
    assert base == pytest.approx(0.993287, abs=2e-6)

"""load_cce_let_table's cce_interp accepts arrays; process_mc_ensemble uses one call.

Bit-identical to the per-event loop (np.interp is elementwise); a scalar still gives a
float, and a plain scalar-only callable still works through the loop.
"""

from __future__ import annotations

import numpy as np
import pandas as pd

from etna.core.mc_coupling import process_mc_ensemble
from etna.core.single_particle import load_cce_let_table


def _events(n=2000, seed=1):
    rng = np.random.default_rng(seed)
    return pd.DataFrame(
        {
            "event_id": np.repeat(np.arange(n), 2),
            "edep_keV": rng.lognormal(3.0, 1.5, 2 * n),
        }
    )


def test_array_lookup_equals_scalar_lookup():
    interp, _ = load_cce_let_table("data/cce_let_table_100um.json")
    lets = np.geomspace(0.05, 5000.0, 997)
    vec = interp(lets)
    assert isinstance(interp(15.0), float)
    assert np.array_equal(vec, np.array([interp(v) for v in lets]))


def test_process_mc_ensemble_vectorized_equals_loop():
    interp, _ = load_cce_let_table("data/cce_let_table_100um.json")
    ev = _events()

    def scalar_only(let):  # no `vectorized` attribute: the per-event loop
        return interp(float(let))

    scalar_only.let_min, scalar_only.let_max = interp.let_min, interp.let_max
    a = process_mc_ensemble(ev, interp)
    b = process_mc_ensemble(ev, scalar_only)
    assert np.array_equal(a["event_CCE"], b["event_CCE"])
    assert np.array_equal(a["event_collected_keV"], b["event_collected_keV"])
    assert a["n_clamped"] == b["n_clamped"]

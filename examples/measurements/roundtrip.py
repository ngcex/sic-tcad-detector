"""Round trip of the SYNTHETIC example measurements (examples/measurements/data/).

Loads each file with the public loaders, compares it with the truth device it was
generated from and rebuilds that device from the SIMS profile. A healthy install
prints residuals at the solver / resampling level; a broken loader, unit or sign
convention shows up as an order-of-magnitude residual. Seconds on the batched engine.

    uv run python examples/measurements/roundtrip.py
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

from etna import measurements as M
from etna.api.device import measured_cv_tabulated_config

DATA = Path(__file__).resolve().parent / "data"


def run() -> dict:
    truth = measured_cv_tabulated_config()
    out = {}

    cv = M.load_cv_csv(DATA / "cv_synthetic.csv")
    cmp = M.compare_cv(cv, truth)
    out["cv"] = {
        "n": len(cv.V),
        "V_min": float(cv.V.min()),
        "rms_W_um": cmp.rms_W_um,
        "rms_C_rel": cmp.rms_C_rel,
    }

    prof = M.load_sims_csv(DATA / "sims_synthetic.csv")
    table = M.resample_profile(prof, truth, method="interp")
    pvm = M.profile_vs_model(table, truth)
    rebuilt = table.to_device_config(truth)
    cmp2 = M.compare_cv(cv, rebuilt)
    out["sims"] = {
        "junction_um": table.junction_depth_um,
        "junction_source": table.junction_source,
        "rms_log10_vs_truth": pvm["rms_log10"],
        "max_abs_log10_vs_truth": pvm["max_abs_log10"],
        "cv_rms_W_um_rebuilt": cmp2.rms_W_um,
        "N_A_rebuilt": rebuilt.N_A,
    }

    it = M.load_it_csv(DATA / "it_synthetic.csv")
    an = M.analyze_it(it, truth)
    out["it"] = {
        f"{c.bias_V:g} V": {
            "E_a_data_eV": c.data_fit.E_a_eV,
            "E_a_model_eV": c.model_fit.E_a_eV,
            "delta_eV": c.delta_eV,
        }
        for c in an.comparisons
    }
    return out


def main():
    res = run()
    c = res["cv"]
    print(
        f"C-V : {c['n']} points to {c['V_min']:g} V; vs truth rms W "
        f"{c['rms_W_um']:.2e} um, rms C {c['rms_C_rel']:.2e}"
    )
    s = res["sims"]
    print(
        f"SIMS: junction {s['junction_um']:g} um ({s['junction_source']}); donor "
        f"table vs truth rms {s['rms_log10_vs_truth']:.3f} dex (max "
        f"{s['max_abs_log10_vs_truth']:.3f}); C-V of the rebuilt device vs the file: "
        f"rms W {s['cv_rms_W_um_rebuilt']:.3g} um"
    )
    for b, d in res["it"].items():
        print(
            f"I(T): {b}: E_a data {d['E_a_data_eV']:.4f} eV, model "
            f"{d['E_a_model_eV']:.4f} eV, delta {d['delta_eV']:+.1e} eV"
        )
    return res


if __name__ == "__main__":
    np.set_printoptions(precision=4)
    main()

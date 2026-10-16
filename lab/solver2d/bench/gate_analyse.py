"""Post-analysis of the gate CSVs (notes/gate_devsim.csv, notes/gate_auto.csv)."""

import numpy as np
import pandas as pd

KEY = ["half_width_um", "epi_thickness_cm", "N_D_bulk", "V_bias"]
d = pd.read_csv("lab/solver2d/notes/gate_devsim.csv")
a = pd.read_csv("lab/solver2d/notes/gate_auto.csv")
m = d.merge(a, on=KEY, suffixes=("_d", "_a"))
ok = m["center_cce_d"].notna() & m["center_cce_a"].notna()
for c in ("center_cce", "edge_cce", "cce_std", "edge_center_ratio"):
    diff = (m[f"{c}_a"] - m[f"{c}_d"])[ok]
    print(f"{c}: max|d| {diff.abs().max():.2e}  signed mean {diff.mean():+.2e}")
print("strict ranking identical:", d[KEY].equals(a[KEY]))


def tie_order(df):
    g = df.copy()
    g["r5"] = g["edge_center_ratio"].round(5)
    g = g.sort_values(KEY).reset_index(drop=True)  # grid order as the tiebreak
    g = g.sort_values(
        ["is_valid", "r5"], ascending=[False, False], na_position="last", kind="stable"
    )
    return g[KEY].reset_index(drop=True)


print(
    "tie-aware ranking (ratio to 5 decimals) identical:",
    tie_order(d).equals(tie_order(a)),
)
print(
    "valid-group ratio spread devsim",
    np.ptp(d.loc[d.is_valid, "edge_center_ratio"]),
    "auto",
    np.ptp(a.loc[a.is_valid, "edge_center_ratio"]),
)
nan = m[m["center_cce_d"].isna() | m["center_cce_a"].isna()]
print(
    nan[
        KEY
        + ["solver_a", "t_config_s_d", "t_config_s_a", "center_cce_d", "center_cce_a"]
    ].to_string()
)
for s in ("d", "a"):
    t = m[f"t_config_s_{s}"]
    print(
        s,
        "per-config time: median",
        round(t.median(), 3),
        "max",
        round(t.max(), 3),
        "sum",
        round(t.sum(), 2),
    )
print(m.groupby("half_width_um")[["t_config_s_d", "t_config_s_a"]].median())

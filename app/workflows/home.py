"""Home landing page: orientation line + current device config summary.

All Streamlit calls live inside render() — no module-level side effects —
so this module can be imported for st.Page registration (app/main.py) and
exercised headlessly via AppTest.from_function(render) in tests.
"""

from __future__ import annotations

import streamlit as st

from app.components.page_utils import require_config


PAGES_TEXT = """
| Page | Computes | Status |
|---|---|---|
| Simulate → Figures | any observable below + publication figures, CSV, animation | as the observable |
| Run All | every simulation for the sidebar device, one table | as each section |
| C-V Analysis | small-signal C(V), 1/C², depletion width | default doping **fitted** to the measured C-V (fit quality, not a prediction) |
| Charge Collection | CCE vs reverse bias, Am-241 5.486 MeV alpha (NIST ASTAR Bragg profile, 18.43 µm range) | predicted, not validated against measurement |
| Field Map | E(x), potential, doping; 2D maps | predicted |
| Radiation Damage | CCE vs proton fluence (SR-NIEL κ(E), carrier removal or charged deep traps) | absolute fluence scale unvalidated (SiC neutron normalisation open) |
| Dark Current | total leakage vs temperature | **calibrated** (N_t fitted to < 18 pA at −60 V): a budget estimate |
| Microdosimetry | lineal-energy spectrum from uploaded MC events (Cauchy mean chord) | geometry-only, y in SiC |
| Batch Sweep | one device parameter swept, curves overlaid | as the observable |
| Export | spec sheet, cross-section, mesh profile, REFERENCE-NOT-FAB layout | layout is not a fabrication layout |
| Measurements | your C-V / SIMS / I(T) files vs the model, C-V inversion | data-driven |
"""


def render() -> None:
    st.title("ETNA — SiC TCAD Simulator")
    st.write(
        "4H-SiC p⁺/n⁻/n⁺ detector simulator (drift-diffusion). Configure a device in "
        "the sidebar, then choose a workflow from the navigation; the sidebar device "
        "is shared by every page."
    )
    st.caption(
        "Solvers: a fast batched 1D engine (C-V, CCE, field, damage, dark current, "
        "FLASH, temperature) and a 2D engine (field maps) inside the parameter box where "
        "they were cross-checked against devsim, devsim (reference finite-volume solver) "
        "otherwise. Cross-checked means the two solvers agree, not that the model is "
        "validated against measurement. Every result shows the solver used."
    )
    st.markdown(PAGES_TEXT)

    cfg = require_config()

    st.subheader("Current device (DeviceConfig)")
    st.caption(
        "Lengths in µm, doping in cm⁻³, T in K, area in cm², N_t in cm⁻³ s⁻¹, E_t in eV "
        "below E_c."
    )
    st.json({k: getattr(cfg, k) for k in cfg.__dataclass_fields__})

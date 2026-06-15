"""Phase 26 — freeze v3.0 reference outputs to tests/baselines/v3_frozen.json.

Computes immutable v3.0 reference values directly via the same Python API the
20 validated notebooks use (NOT by re-running notebooks) and writes them to
``tests/baselines/v3_frozen.json``:

  - ``cce_center_100um`` / ``cce_center_300um`` : whole-device CCE (uniform epi
    generation) for the 100 um and 300 um sensitive volumes at a low reverse
    bias (V_bias = 10 V, well below the -15 V v3.0 convergence ceiling).
  - ``cv_1d`` : 1D depletion width W(V) at 0 / -5 / -10 V (the v3.0 1D-calibrated
    targets the 2D solver must reproduce at the centerline in Plan 03).
  - ``notebook_list`` : the canonical 20 scientific-deliverable notebooks.
  - ``metadata`` : freeze timestamp, git SHA (or WORKING_TREE), tolerance
    (0.1 %), devsim version, v3.0 default doping, purpose.

Plan 04 asserts that the calibrated 2D defaults reproduce these within
``metadata.tolerance_rel`` at biases <= -10 V — any src/ change that drifts the
values by > 0.1 % fails CI.

DOES NOT modify src/. DOES NOT re-run any notebook.

Run:
    cd /Users/ngcex/projects/physics/petringa && uv run python scripts/freeze_v3_baselines.py

NOTE on devsim global state: ``devsim.reset_devsim()`` clears the direct-solver
parameter AND the UMFPACK solver callback; this script saves and restores both
between device builds (mirrors src/optimization.py), otherwise the second build
fails with ``Solver "custom" specified, but "solver_callback" not set``.
"""

import datetime
import json
import pathlib
import subprocess
import sys
import uuid

# Make `import src.*` resolve when run as `uv run python scripts/...`.
_PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import devsim
import numpy as np

from src.charge_collection_2d import (
    create_2d_dd_device,
    compute_cce_2d,
    _robust_dc_solve,
)
from src.charge_collection import add_generation_to_dd
from src.drift_diffusion import create_dd_device
from src.cv_analysis import cv_sweep

# v3.0 default graded doping (src/device2d.py:41-43)
N_D_JUNCTION = 2.9e15
N_D_BULK = 8.5e13
L_TRANSITION = 1e-4

# Canonical 20 scientific-deliverable notebooks (research A3 / D4: the duplicates
# 03_executed.ipynb and 05_parametric_studies.ipynb are NOT deliverables).
NOTEBOOK_LIST = [
    "01_phase1_validation.ipynb",
    "02_electrical_characterization.ipynb",
    "03_charge_collection.ipynb",
    "04_flash_recombination.ipynb",
    "05_dark_current_vs_fluence.ipynb",
    "06_temperature_dependence.ipynb",
    "07_dark_current.ipynb",
    "08_transient_flash.ipynb",
    "09_radiation_damage.ipynb",
    "10_cce_vs_fluence.ipynb",
    "11_dark_current_cv_evolution.ipynb",
    "12_multi_defect_comparison.ipynb",
    "13_parametric_optimization.ipynb",
    "14_validation.ipynb",
    "15_2d_electrostatics_cce.ipynb",
    "16_single_particle_cce.ipynb",
    "17_mc_coupling.ipynb",
    "18_microdosimetric_spectra.ipynb",
    "19_alternative_structures.ipynb",
    "20_feasibility_report.ipynb",
]

OUTPUT_PATH = pathlib.Path("tests/baselines/v3_frozen.json")


def _reset_devsim_safe():
    """Reset devsim and restore direct_solver + solver_callback (PITFALLS P03)."""
    try:
        saved_solver = devsim.get_parameter(name="direct_solver")
    except Exception:
        saved_solver = "custom"
    try:
        saved_callback = devsim.get_parameter(name="solver_callback")
    except Exception:
        saved_callback = None
    devsim.reset_devsim()
    try:
        devsim.set_parameter(name="direct_solver", value=saved_solver)
        if saved_callback is not None:
            devsim.set_parameter(name="solver_callback", value=saved_callback)
    except Exception:
        pass


def compute_cce_center(half_width_um):
    """Whole-device CCE with uniform epi generation at V_bias = 10 V.

    Uses the proven test recipe (tests/test_charge_collection_2d.py:120-126):
    build the 2D DD device at a low reverse bias that v3.0 converges at, inject a
    uniform generation rate across the epi (y >= junction_pos), solve, and read
    the CCE scalar from compute_cce_2d.
    """
    _reset_devsim_safe()
    name = f"freeze2d_{uuid.uuid4().hex[:8]}"
    try:
        dev = create_2d_dd_device(
            device_name=name,
            half_width_um=half_width_um,
            V_bias=10.0,  # below the -15 V v3.0 ceiling; v3.0 converges here
            doping_profile="graded",
            N_D_junction=N_D_JUNCTION,
            N_D_bulk=N_D_BULK,
            L_transition=L_TRANSITION,
        )
        device = dev["device_name"]
        region = dev["region_name"]
        junction_pos = dev["substrate_thickness_cm"]
        y_nodes = np.array(
            devsim.get_node_model_values(device=device, region=region, name="y")
        )
        gen_rate = 1e18  # cm^-3 s^-1, uniform across epi
        gen = np.where(y_nodes >= junction_pos, gen_rate, 0.0)
        add_generation_to_dd(dev, gen)
        _robust_dc_solve()
        cce = compute_cce_2d(dev, gen)
        return float(cce)
    finally:
        try:
            devsim.delete_device(device=name)
        except Exception:
            pass


def compute_cv_1d():
    """1D depletion width W(V) at 0 / -5 / -10 V via cv_sweep (graded profile)."""
    _reset_devsim_safe()
    name = f"freeze1d_{uuid.uuid4().hex[:8]}"
    try:
        dev = create_dd_device(
            device_name=name,
            doping_profile="graded",
            N_D_junction=N_D_JUNCTION,
            N_D_bulk=N_D_BULK,
            L_transition=L_TRANSITION,
        )
        result = cv_sweep(dev, V_range=[0.0, -5.0, -10.0])
        voltages = list(result["voltages"])
        widths = list(result["depletion_widths"])
        cv = {}
        for v, w in zip(voltages, widths):
            cv[f"{float(v)}"] = float(w)
        # Guarantee the three required keys exist (string form matching V_range).
        for key in ("0.0", "-5.0", "-10.0"):
            cv.setdefault(key, None)
        return cv
    finally:
        try:
            devsim.delete_device(device=name)
        except Exception:
            pass


def build_metadata():
    porcelain = (
        subprocess.check_output(["git", "status", "--porcelain"]).decode().strip()
    )
    if porcelain == "":
        git_sha = subprocess.check_output(["git", "rev-parse", "HEAD"]).decode().strip()
    else:
        git_sha = "WORKING_TREE"
    return {
        "frozen_at_utc": datetime.datetime.now(datetime.timezone.utc).isoformat(),
        "git_commit_sha": git_sha,
        "tolerance_rel": 0.001,
        "devsim_version": getattr(devsim, "__version__", "unknown"),
        "v3_default_doping": {
            "N_D_junction": N_D_JUNCTION,
            "N_D_bulk": N_D_BULK,
            "L_transition": L_TRANSITION,
        },
        "purpose": (
            "Phase 26 regression baseline — Plan 04 asserts that calibrated 2D "
            "defaults reproduce these within tolerance_rel at biases <= -10 V"
        ),
    }


def main():
    # Ensure all canonical notebooks exist before freezing.
    missing = [
        nb for nb in NOTEBOOK_LIST if not pathlib.Path(f"notebooks/{nb}").exists()
    ]
    if missing:
        raise SystemExit(f"Missing canonical notebooks: {missing}")

    pathlib.Path("tests/baselines").mkdir(parents=True, exist_ok=True)

    print("[freeze] computing CCE (100 um SV)...")
    cce_100 = compute_cce_center(half_width_um=50.0)
    print(f"[freeze]   cce_center_100um = {cce_100:.6f}")

    print("[freeze] computing CCE (300 um SV)...")
    cce_300 = compute_cce_center(half_width_um=150.0)
    print(f"[freeze]   cce_center_300um = {cce_300:.6f}")

    print("[freeze] computing 1D C-V at 0 / -5 / -10 V...")
    cv_1d = compute_cv_1d()
    print(f"[freeze]   cv_1d = {cv_1d}")

    payload = {
        "metadata": build_metadata(),
        "cce_center_100um": cce_100,
        "cce_center_300um": cce_300,
        "cv_1d": cv_1d,
        "notebook_list": NOTEBOOK_LIST,
    }

    with open(OUTPUT_PATH, "w") as f:
        json.dump(payload, f, indent=2, sort_keys=True)
    print(f"[freeze] wrote {OUTPUT_PATH}")
    return 0


if __name__ == "__main__":
    sys.exit(main())

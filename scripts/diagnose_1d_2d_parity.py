"""Phase 26 — 1D vs 2D reverse-bias convergence diagnostic.

Standalone, runnable diagnostic that ramps a 1D and two 2D devices (100 um and
300 um sensitive volumes) with the *identical* v3.0 graded doping profile
``{N_D_junction=2.9e15, N_D_bulk=8.5e13, L_transition=1e-4}`` from 0 V toward
-50 V reverse bias and records the last converged cathode voltage for each.

It then classifies the root cause of the known ">-15 V solver divergence" as one
of three hypotheses and writes the result to
``.planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md`` with a
machine-readable ``hypothesis:`` YAML key (one of H1/H2/H3):

  - H1 (profile)        : the graded profile itself is marginal (1D also fails
                          before -30 V).
  - H2 (mesh/BC)        : 1D reaches <= -30 V but 2D fails 2D-specifically before
                          -15 V (mesh / boundary-condition differential).
  - H3 (both)           : 1D fine, 2D reaches between -15 V and -30 V — both the
                          profile margin and the mesh-BC differential contribute.

This artifact is consumed by Plan 03 (26-03-PLAN.md Task 1) to gate the
calibration approach: jumping straight to a 2D Nelder-Mead refit only makes sense
if H1/H3 hold; if H2 dominates, the fix is mesh/BC, not the profile.

DIAGNOSTIC ONLY — this script does NOT modify any file under ``src/``.

Run:
    cd /Users/ngcex/projects/physics/petringa && uv run python scripts/diagnose_1d_2d_parity.py

NOTE on imports: ``create_dd_device`` is exported from ``src.drift_diffusion``
(NOT ``src.device`` — the plan's import path was corrected here; ``src.device``
only exposes the lower-level ``create_sic_device``).
"""

import datetime
import pathlib
import sys
import uuid

# Allow `uv run python scripts/diagnose_1d_2d_parity.py` to find the `src`
# package: when run as a script, sys.path[0] is the scripts/ dir, not the
# project root. Prepend the project root so `import src.*` resolves without
# requiring PYTHONPATH to be set externally.
_PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import devsim
import numpy as np

from src.drift_diffusion import create_dd_device
from src.charge_collection_2d import create_2d_dd_device
from src.poisson import extract_depletion_width_numerical


def _reset_devsim_safe():
    """Reset devsim and restore the direct-solver settings.

    ``devsim.reset_devsim()`` clears BOTH the ``direct_solver`` parameter (to the
    sentinel ``"unknown"``) and the ``solver_callback`` that registers the
    UMFPACK ``custom`` solver. Without restoring both, the next ``devsim.solve``
    raises either ``Unrecognized "direct_solver" parameter value "unknown"`` or
    ``Solver "custom" specified, but "solver_callback" not set``. Save and
    restore both — this mirrors the proven pattern in ``src/optimization.py``
    (lines 112-120, 199-205). (Blocking-issue fix; see PITFALLS P03 on devsim
    global-state hygiene.)
    """
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


# --- Constants: copied verbatim from src/device2d.py:41-43 (v3.0 defaults) ---
N_D_JUNCTION = 2.9e15  # cm^-3, near junction
N_D_BULK = 8.5e13  # cm^-3, deep epi bulk
L_TRANSITION = 1e-4  # cm, characteristic decay length
EPI_THICKNESS_CM = 10e-4
SUBSTRATE_THICKNESS_CM = 1e-4
N_A = 1e19
T = 300

# Reverse-bias targets (conventional negative sign). The script ramps the cathode
# in 0.5 V steps and records the last converged reverse-bias voltage.
TARGET_VOLTAGES = [-5.0, -10.0, -15.0, -20.0, -30.0, -40.0, -50.0]

OUTPUT_PATH = pathlib.Path(
    ".planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md"
)


def _unique(prefix):
    return f"{prefix}_{uuid.uuid4().hex[:8]}"


def ramp_to_max(create_fn, device_label, **kwargs):
    """Build a device via ``create_fn`` then ramp the cathode toward -50 V.

    Ramps in 0.5 V steps using the proven cv_sweep fallback tolerances
    (absolute_error=1e10, relative_error=1e-10, max_iter=40; on devsim.error
    retry with 1e12, 1e-8, 100). Tracks the last cathode voltage at which
    ``devsim.solve`` succeeded.

    Returns
    -------
    (V_max_reverse_bias, failure_message)
        ``V_max_reverse_bias`` is reported as a NEGATIVE number (conventional
        reverse bias). On total failure at the first step returns
        ``(0.0, "did not converge below 0 V")``. The device is always deleted in
        a ``finally`` block (PITFALLS P20 — unique names, clean teardown).
    """
    name = kwargs.get("device_name")
    if name is None:
        name = _unique("diag")
        kwargs["device_name"] = name

    last_converged_reverse = 0.0
    failure_message = "reached -50 V without failure"

    try:
        # Build the device with the cathode at 0 V (V_bias=0). The 2D wrapper
        # ramps to V_bias internally; passing 0.0 means "build at equilibrium".
        device_info = create_fn(**kwargs)

        # Cathode bias name. For p+/n- diode, positive cathode V = reverse bias.
        try:
            from src import cv_analysis as _cv  # for GetContactBiasName

            bias_name = _cv.simple_physics.GetContactBiasName("cathode")
        except Exception:
            bias_name = "cathode_bias"

        device = device_info["device_name"]

        V_cathode = 0.0
        step = 0.5
        target_cathode = 50.0  # corresponds to -50 V reverse bias

        while V_cathode < target_cathode - 1e-9:
            V_cathode = round(V_cathode + step, 6)
            reverse_bias = -V_cathode
            print(f"[diag] {device_label} ramping cathode to {reverse_bias:.2f} V")
            devsim.set_parameter(device=device, name=bias_name, value=V_cathode)
            try:
                devsim.solve(
                    type="dc",
                    absolute_error=1e10,
                    relative_error=1e-10,
                    maximum_iterations=40,
                )
            except devsim.error:
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e12,
                        relative_error=1e-8,
                        maximum_iterations=100,
                    )
                except devsim.error as e:
                    failure_message = (
                        f"diverged at {reverse_bias:.2f} V: {str(e)[:160]}"
                    )
                    break
            last_converged_reverse = reverse_bias

        if last_converged_reverse == 0.0:
            return 0.0, "did not converge below 0 V"
        return last_converged_reverse, failure_message
    except Exception as e:  # build-time or unexpected failure
        return 0.0, f"build/ramp error: {str(e)[:160]}"
    finally:
        try:
            devsim.delete_device(device=name)
        except Exception:
            pass


def _W_at_voltages_1d(voltages):
    """Build a fresh 1D device, ramp to each voltage, return {V: W_cm}.

    Each voltage gets an isolated build (devsim.reset_devsim before each) so the
    extraction is unaffected by prior ramp state. Returns None for any voltage
    that fails to converge.
    """
    out = {}
    from src import cv_analysis as _cv

    bias_name = _cv.simple_physics.GetContactBiasName("cathode")
    for v in voltages:
        _reset_devsim_safe()
        name = _unique("diag1dW")
        try:
            dev = create_dd_device(
                device_name=name,
                doping_profile="graded",
                N_D_junction=N_D_JUNCTION,
                N_D_bulk=N_D_BULK,
                L_transition=L_TRANSITION,
                epi_thickness_cm=EPI_THICKNESS_CM,
                substrate_thickness_cm=SUBSTRATE_THICKNESS_CM,
                N_A=N_A,
                T=T,
            )
            device = dev["device_name"]
            V_cathode = 0.0
            target = -v  # conventional reverse bias -> positive cathode
            converged = True
            while V_cathode < target - 1e-9:
                V_cathode = round(V_cathode + 0.5, 6)
                devsim.set_parameter(device=device, name=bias_name, value=V_cathode)
                try:
                    devsim.solve(
                        type="dc",
                        absolute_error=1e10,
                        relative_error=1e-10,
                        maximum_iterations=40,
                    )
                except devsim.error:
                    try:
                        devsim.solve(
                            type="dc",
                            absolute_error=1e12,
                            relative_error=1e-8,
                            maximum_iterations=100,
                        )
                    except devsim.error:
                        converged = False
                        break
            out[v] = (
                float(extract_depletion_width_numerical(dev)) if converged else None
            )
        except Exception:
            out[v] = None
        finally:
            try:
                devsim.delete_device(device=name)
            except Exception:
                pass
    return out


def run_1d():
    """Ramp a 1D graded device and extract W at -10 V and -30 V."""
    _reset_devsim_safe()
    V_max, fail = ramp_to_max(
        create_dd_device,
        "1D",
        device_name=_unique("diag1d"),
        doping_profile="graded",
        N_D_junction=N_D_JUNCTION,
        N_D_bulk=N_D_BULK,
        L_transition=L_TRANSITION,
        epi_thickness_cm=EPI_THICKNESS_CM,
        substrate_thickness_cm=SUBSTRATE_THICKNESS_CM,
        N_A=N_A,
        T=T,
    )
    # W extraction (1D only — extract_depletion_width_numerical reads x as depth)
    W = _W_at_voltages_1d([-10.0, -30.0])
    return {
        "V_max": V_max,
        "failure_message": fail,
        "W_at_minus_10V_cm": W.get(-10.0),
        "W_at_minus_30V_cm": W.get(-30.0),
    }


def run_2d(half_width_um):
    """Ramp a 2D graded device of the given half-width toward -50 V.

    W extraction is intentionally NOT performed: extract_depletion_width_numerical
    is 1D-only (reads x as depth). The 2D-aware center-column extractor is Plan 02
    work, so W is reported as None for 2D here.
    """
    _reset_devsim_safe()
    V_max, fail = ramp_to_max(
        create_2d_dd_device,
        f"2D {int(half_width_um * 2)}um SV",
        device_name=_unique("diag2d"),
        half_width_um=half_width_um,
        V_bias=0.0,  # build at equilibrium; ramp_to_max does the sweep
        doping_profile="graded",
        N_D_junction=N_D_JUNCTION,
        N_D_bulk=N_D_BULK,
        L_transition=L_TRANSITION,
        epi_thickness_cm=EPI_THICKNESS_CM,
        substrate_thickness_cm=SUBSTRATE_THICKNESS_CM,
        N_A=N_A,
        T=T,
    )
    return {
        "V_max": V_max,
        "failure_message": fail,
        "W_at_minus_10V_cm": None,
        "W_at_minus_30V_cm": None,
    }


def classify(d1, d2_100, d2_300):
    """Classify root cause as H1 / H2 / H3 per deterministic rules.

    Voltages are negative for reverse bias; "more negative = deeper bias".
    """
    d2_min = min(d2_100["V_max"], d2_300["V_max"])  # the *worse* (less negative) 2D
    # 1D itself fails before -30 V -> profile is marginal
    if d1["V_max"] > -30.0:
        return "H1"
    # 1D OK to <= -30 V but 2D fails 2D-specifically before -15 V -> mesh/BC
    if d1["V_max"] <= -30.0 and d2_min > -15.0:
        return "H2"
    # 1D OK; 2D reaches between -15 V and -30 V -> both contribute
    if d1["V_max"] <= -30.0 and d2_min <= -15.0 and d2_min > -30.0:
        return "H3"
    # 1D fine, 2D reaches -30 V or beyond -> profile already adequate
    return "H1"


_EXPLANATIONS = {
    "H1": (
        "Both the 1D and/or 2D devices fail to reach -30 V with the v3.0 "
        "`{2.9e15, 8.5e13, 1e-4}` profile, OR both reach -30 V comfortably. In "
        "either reading the *profile* is the dominant lever: if it fails, the "
        "profile is marginal and a Nelder-Mead refit (Plan 03) is justified; if "
        "both converge deep, the profile is already adequate and calibration is "
        "fine-tuning only. Either way the lever is the doping profile, not the "
        "mesh/BC."
    ),
    "H2": (
        "The 1D device with the identical graded profile converges to at least "
        "-30 V, but the 2D device(s) diverge before -15 V. Because the profile is "
        "demonstrably fine in 1D, the failure is a 2D-specific mesh / "
        "boundary-condition differential (symmetry BC at x=0, air-buffer "
        "interaction near the cathode, or lateral x-mesh coupling). A pure "
        "profile refit will NOT fix this — Plan 03 must address mesh/BC."
    ),
    "H3": (
        "The 1D device converges to at least -30 V, but the 2D device(s) reach "
        "only between -15 V and -30 V. Both the profile margin and a 2D-specific "
        "mesh/BC differential contribute. Plan 03's calibration should combine a "
        "profile refit with attention to the 2D mesh/BC envelope."
    ),
}


def _fmt(v):
    return "n/a" if v is None else f"{v:.3e}"


def main():
    print("[diag] Phase 26 1D-vs-2D reverse-bias convergence diagnostic")
    print("[diag] profile = {2.9e15, 8.5e13, 1e-4} (v3.0 defaults)")

    try:
        d1 = run_1d()
    except devsim.error as e:
        d1 = {
            "V_max": 0.0,
            "failure_message": str(e)[:160],
            "W_at_minus_10V_cm": None,
            "W_at_minus_30V_cm": None,
        }

    try:
        d2_100 = run_2d(half_width_um=50.0)  # 100 um SV
    except devsim.error as e:
        d2_100 = {
            "V_max": 0.0,
            "failure_message": str(e)[:160],
            "W_at_minus_10V_cm": None,
            "W_at_minus_30V_cm": None,
        }

    try:
        d2_300 = run_2d(half_width_um=150.0)  # 300 um SV
    except devsim.error as e:
        d2_300 = {
            "V_max": 0.0,
            "failure_message": str(e)[:160],
            "W_at_minus_10V_cm": None,
            "W_at_minus_30V_cm": None,
        }

    hypothesis = classify(d1, d2_100, d2_300)
    ts = datetime.datetime.now(datetime.timezone.utc).isoformat()

    body = f"""---
phase: 26
document: DIAGNOSIS
created: {ts}
hypothesis: {hypothesis}
---

# Phase 26 — 1D vs 2D Reverse-Bias Convergence Diagnosis

Profile under test (v3.0 defaults): `N_D_junction=2.9e15`, `N_D_bulk=8.5e13`,
`L_transition=1e-4` cm. Each device was ramped from 0 V toward -50 V in 0.5 V
cathode steps using the cv_sweep adaptive-fallback solver pattern.

## Results

| Device          | V_max_reverse_bias (V) | W(-10 V) (cm) | W(-30 V) (cm) | Failure mode |
|-----------------|------------------------|---------------|---------------|--------------|
| 1D (device.py)  | {d1['V_max']:.2f}             | {_fmt(d1['W_at_minus_10V_cm'])}    | {_fmt(d1['W_at_minus_30V_cm'])}    | {d1['failure_message']}    |
| 2D 100 µm SV    | {d2_100['V_max']:.2f}         | n/a           | n/a           | {d2_100['failure_message']}|
| 2D 300 µm SV    | {d2_300['V_max']:.2f}         | n/a           | n/a           | {d2_300['failure_message']}|

Machine-readable evidence:

- `V_1d_max_converged`: {d1['V_max']:.2f} V
- `V_2d_max_converged_100um`: {d2_100['V_max']:.2f} V
- `V_2d_max_converged_300um`: {d2_300['V_max']:.2f} V

## Hypothesis: {hypothesis}

{_EXPLANATIONS[hypothesis]}

## Consumed by

- .planning/phases/26-graded-doping-2d-calibration/26-03-PLAN.md (Task 1 reads `hypothesis:` and gates calibration)
"""

    OUTPUT_PATH.parent.mkdir(parents=True, exist_ok=True)
    OUTPUT_PATH.write_text(body)
    print(f"[diag] wrote {OUTPUT_PATH} (hypothesis={hypothesis})")
    return 0


if __name__ == "__main__":
    sys.exit(main())

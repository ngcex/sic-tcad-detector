#!/usr/bin/env python
"""Phase 26 — canonical end-to-end runner for the graded-doping 2D calibration.

This is the single entry point that ties together the Phase 26 calibration:

  1. Reads the root-cause hypothesis from
     ``.planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md``
     (produced by Plan 01 Task 1).
  2. Gates on the hypothesis: H1 or H3 -> proceed with the Nelder-Mead refit;
     H2 -> abort (exit code 2) because the root cause is mesh / boundary-condition
     / lateral-mesh, not the doping profile, so a profile refit would be wasted
     effort.
  3. Runs :func:`src.device2d.calibrate_graded_doping_2d` on the 100 um SV.
  4. Writes the optimised parameters + per-voltage W + final cost to
     ``26-CALIBRATION-RESULT.md`` (frozen record).
  5. Patches the three ``_N_D_*_DEFAULT`` constants in ``src/device2d.py``
     in-place (line-anchored regex -> idempotent), tagging them
     ``# calibrated Plan 26-03``.

Run (default maxiter=80; ~25-50 min):
    uv run python scripts/run_calibration_2d.py

Options:
    --maxiter N        Nelder-Mead iteration cap (default 80).
    --half_width_um F  SV half-width in um (default 50.0 -> 100 um SV).
    --skip-patch       Run + write result MD but do NOT edit src/device2d.py.
"""

import argparse
import datetime
import pathlib
import re
import sys

# Allow `uv run python scripts/run_calibration_2d.py` to find the `src`
# package: when run as a script, sys.path[0] is the scripts/ dir, not the
# project root. Prepend the project root so `import src.*` resolves without
# requiring PYTHONPATH to be set externally.
_PROJECT_ROOT = pathlib.Path(__file__).resolve().parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from src.device2d import calibrate_graded_doping_2d  # noqa: E402

DIAGNOSIS_PATH = pathlib.Path(
    ".planning/phases/26-graded-doping-2d-calibration/26-DIAGNOSIS.md"
)
RESULT_PATH = pathlib.Path(
    ".planning/phases/26-graded-doping-2d-calibration/26-CALIBRATION-RESULT.md"
)
DEVICE2D_PATH = pathlib.Path("src/device2d.py")

# v3.0 starting values (for the result-MD comparison columns).
_V30_N_D_JUNCTION = 2.9e15
_V30_N_D_BULK = 8.5e13
_V30_L_TRANSITION = 1.0e-4

# Petringa 1D-twin W targets (cm) — the calibration target at the device center.
_W_TARGETS = {0.0: 1.7e-4, -10.0: 9.5e-4, -30.0: 9.73e-4}


def read_hypothesis():
    """Read the ``hypothesis:`` value (H1/H2/H3) from 26-DIAGNOSIS.md."""
    if not DIAGNOSIS_PATH.exists():
        print(
            f"ERROR: {DIAGNOSIS_PATH} not found. Run Plan 01 "
            "(scripts/diagnose_1d_2d_parity.py) first to produce the diagnosis.",
            file=sys.stderr,
        )
        sys.exit(1)

    text = DIAGNOSIS_PATH.read_text()
    match = re.search(r"^hypothesis:\s*(H1|H2|H3)\s*$", text, flags=re.MULTILINE)
    if match is None:
        print(
            "ERROR: DIAGNOSIS.md exists but contains no parseable "
            "'hypothesis:' line (expected H1, H2, or H3).",
            file=sys.stderr,
        )
        sys.exit(1)

    return match.group(1)


def gate_on_hypothesis(h):
    """Abort on H2 (mesh/BC), proceed on H1/H3 (profile-driven)."""
    if h == "H2":
        print(
            "ABORT: Plan 01 diagnosis classified the root cause as H2 "
            "(2D-specific mesh / boundary-condition differential, NOT the "
            "doping profile).\n"
            "\n"
            "A Nelder-Mead refit of {N_D_junction, N_D_bulk, L_transition} "
            "will NOT fix an H2 failure. The suspect causes are:\n"
            "  1. air-buffer interaction with the near-cathode SCR at high "
            "reverse bias\n"
            "  2. x-symmetry boundary condition behavior in the 2D mesh\n"
            "  3. lateral mesh-line spacing differential between 1D and 2D\n"
            "\n"
            "Investigate 26-DIAGNOSIS.md (and re-run the diagnostic) before "
            "re-running calibration. See 26-RESEARCH.md Open Questions item 1 "
            "for details.",
            file=sys.stderr,
        )
        sys.exit(2)

    print(f"[gate] hypothesis={h} proceeding to calibration")


def run_calibration(maxiter, half_width_um):
    """Run the 2D graded-doping calibration and return the result dict."""
    return calibrate_graded_doping_2d(maxiter=maxiter, half_width_um=half_width_um)


def write_result_md(result):
    """Write the frozen calibration result to RESULT_PATH (YAML + markdown)."""
    created = datetime.datetime.now(datetime.timezone.utc).isoformat()

    # Per-voltage W check rows.
    w_rows = []
    for v in sorted(result["W_simulated"].keys(), reverse=True):
        w_sim = result["W_simulated"][v]
        w_target = _W_TARGETS.get(v)
        if w_target is not None and w_target != 0:
            rel = abs(w_sim - w_target) / w_target
            rel_str = f"{rel:.4f}"
            target_str = f"{w_target:.3e}"
        else:
            rel_str = "n/a"
            target_str = "n/a"
        w_rows.append(f"| {v:.1f} | {w_sim:.3e} | {target_str} | {rel_str} |")
    w_table = "\n".join(w_rows)

    frontmatter = (
        "---\n"
        "phase: 26\n"
        "document: CALIBRATION_RESULT\n"
        f"created: {created}\n"
        f"N_D_junction: {result['N_D_junction']:.6e}\n"
        f"N_D_bulk: {result['N_D_bulk']:.6e}\n"
        f"L_transition: {result['L_transition']:.6e}\n"
        f"final_cost: {result['final_cost']:.6e}\n"
        f"nit: {result['nit']}\n"
        f"success: {result['success']}\n"
        f"converged_at_minus_50V: {result['converged_at_convergence_target']}\n"
        "---\n"
    )

    body = (
        "\n# Phase 26 — Graded Doping 2D Calibration Result\n"
        "\n"
        "## Optimised parameters\n"
        "\n"
        "| Parameter | Calibrated value | v3.0 starting value |\n"
        "| --------- | ---------------- | ------------------- |\n"
        f"| N_D_junction | {result['N_D_junction']:.3e} | {_V30_N_D_JUNCTION:.1e} |\n"
        f"| N_D_bulk | {result['N_D_bulk']:.3e} | {_V30_N_D_BULK:.1e} |\n"
        f"| L_transition cm | {result['L_transition']:.3e} | {_V30_L_TRANSITION:.1e} |\n"
        "\n"
        "## Per-voltage W check (2D center column)\n"
        "\n"
        "| V (V) | W_sim (cm) | W_target (cm) | rel error |\n"
        "| ----- | ---------- | ------------- | --------- |\n"
        f"{w_table}\n"
        "\n"
        "## Convergence\n"
        "\n"
        f"- Nelder-Mead iterations: {result['nit']}\n"
        f"- final_cost: {result['final_cost']:.6e}\n"
        "- convergence target (-50 V) reached: "
        f"{result['converged_at_convergence_target']}\n"
        "\n"
        "## Consumed by\n"
        "\n"
        "- `src/device2d.py` `_N_D_JUNCTION_DEFAULT` / `_N_D_BULK_DEFAULT` / "
        "`_L_TRANSITION_DEFAULT` constants (patched by this script).\n"
        "- `tests/test_device2d.py::TestReverseBiasConvergence` + "
        "`TestCalibrationCV` (Plan 03 Task 3).\n"
        "- `tests/test_v3_baseline_regression.py` (Plan 04 — v3.0 low-bias "
        "regression guard).\n"
    )

    RESULT_PATH.write_text(frontmatter + body)
    print(f"[result] wrote {RESULT_PATH}")


def patch_device2d_defaults(result):
    """Patch the three _N_D_*_DEFAULT constants in src/device2d.py in-place."""
    text = DEVICE2D_PATH.read_text()

    replacements = [
        (
            r"^_N_D_JUNCTION_DEFAULT\s*=\s*[0-9.eE+\-]+.*$",
            f"_N_D_JUNCTION_DEFAULT = {result['N_D_junction']:.6e}  "
            "# cm^-3, calibrated Plan 26-03",
            "_N_D_JUNCTION_DEFAULT",
        ),
        (
            r"^_N_D_BULK_DEFAULT\s*=\s*[0-9.eE+\-]+.*$",
            f"_N_D_BULK_DEFAULT = {result['N_D_bulk']:.6e}  "
            "# cm^-3, calibrated Plan 26-03",
            "_N_D_BULK_DEFAULT",
        ),
        (
            r"^_L_TRANSITION_DEFAULT\s*=\s*[0-9.eE+\-]+.*$",
            f"_L_TRANSITION_DEFAULT = {result['L_transition']:.6e}  "
            "# cm, calibrated Plan 26-03",
            "_L_TRANSITION_DEFAULT",
        ),
    ]

    for pattern, replacement, name in replacements:
        text, count = re.subn(pattern, replacement, text, count=1, flags=re.MULTILINE)
        if count != 1:
            raise SystemExit(
                f"ERROR: failed to patch constant {name} in {DEVICE2D_PATH} "
                f"(matched {count} times, expected 1)."
            )

    DEVICE2D_PATH.write_text(text)
    print(f"[patch] updated {DEVICE2D_PATH} constants with calibrated values")


def main():
    parser = argparse.ArgumentParser(
        description="Run the Phase 26 graded-doping 2D calibration end-to-end."
    )
    parser.add_argument("--maxiter", type=int, default=80)
    parser.add_argument("--half_width_um", type=float, default=50.0)
    parser.add_argument("--skip-patch", action="store_true", default=False)
    args = parser.parse_args()

    h = read_hypothesis()
    gate_on_hypothesis(h)

    result = run_calibration(args.maxiter, args.half_width_um)
    write_result_md(result)

    if not args.skip_patch:
        patch_device2d_defaults(result)

    print(
        "[done] calibration complete: "
        f"N_D_junction={result['N_D_junction']:.3e}, "
        f"N_D_bulk={result['N_D_bulk']:.3e}, "
        f"L_transition={result['L_transition']:.3e} cm, "
        f"final_cost={result['final_cost']:.6e}, "
        f"converged_at_-50V={result['converged_at_convergence_target']}"
    )


if __name__ == "__main__":
    main()

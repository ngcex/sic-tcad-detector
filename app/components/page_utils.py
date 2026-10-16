"""Shared Streamlit helpers for the workflow pages.

Error reporting, stale-result warning, solver/validated-domain badge and the
full-depletion / truncation notes. Text builders (``solver_badge_text``,
``reach_through_text``, ``truncation_text``) are pure so they can be unit
tested; the ``show_*`` functions render them.
"""

from __future__ import annotations

import traceback

import numpy as np
import streamlit as st

SOLVER_LABELS = {
    "batched": "batched 1D engine (numba)",
    "engine2d": "2D engine (numba)",
    "devsim": "devsim (reference finite-volume solver)",
}

STALE_RESULT_TEXT = (
    "These results are for a previous configuration. The sidebar settings have "
    "changed since this run; press **Run simulation** to update."
)


def show_error(exc: BaseException, hint: str = "") -> None:
    """Friendly st.error plus an expander with the technical details."""
    name = type(exc).__name__
    if isinstance(exc, ValueError):
        head = f"Invalid input: {exc}"
    elif isinstance(exc, NotImplementedError):
        head = f"Not supported: {exc}"
    elif isinstance(exc, RuntimeError) or name == "error":  # devsim.error
        head = "The solver did not converge for this device / bias. " + (
            hint or "Try a smaller bias range or different device parameters."
        )
    else:
        head = f"The simulation failed unexpectedly ({name})."
        if hint:
            head += " " + hint
    st.error(head)
    with st.expander("Error details"):
        st.code(
            "".join(traceback.format_exception(type(exc), exc, exc.__traceback__))
            or f"{name}: {exc}"
        )


def show_stale_warning(result_config, cfg) -> None:
    """Warn when the cached result was computed for a different sidebar config."""
    if result_config is not None and cfg is not None and result_config != cfg:
        st.warning(STALE_RESULT_TEXT)


def solver_badge_text(metadata: dict) -> str:
    """'Solver: ... · inside/outside the validated domain' (missing keys tolerated)."""
    solver = metadata.get("solver")
    if solver is None:
        return ""
    label = SOLVER_LABELS.get(solver, solver)
    parts = [f"Solver: {label}"]
    if metadata.get("fallback_reason"):
        parts[0] += f" (fallback: {metadata['fallback_reason']})"
    if "validated_domain" in metadata:
        if metadata["validated_domain"]:
            parts.append("device inside the validated domain")
        else:
            note = metadata.get("domain_note")
            parts.append(
                "device OUTSIDE the validated domain"
                + (f" ({note})" if note else "")
                + "; no engine/devsim cross-check"
            )
    return " · ".join(parts)


def reach_through_text(result) -> str:
    """Full-depletion note (informational) for C-V / CCE results."""
    md = result.metadata
    if "reach_through_V" not in md:
        return ""
    v_rt = md.get("reach_through_V")
    if v_rt is None:
        if md.get("reach_through_known", True):
            return "The epi layer is not fully depleted within this bias range."
        return (
            "Full-depletion voltage not computed (device outside the validated "
            "domain)."
        )
    return (
        f"Epi fully depleted from |V| ≈ {v_rt:.1f} V (dotted line). Beyond it the "
        "capacitance saturates and the result depends on the cathode model "
        "(ohmic contact directly on the n⁻ epi, no n⁺ substrate)."
    )


def cross_check_text(result) -> str:
    """Flag engine points deeper than the engine-vs-devsim parity checks reach."""
    v_cc = result.metadata.get("cross_checked_to_V")
    if v_cc is None:
        return ""
    x = np.asarray(result.x, dtype=float)
    beyond = int(np.sum(np.abs(x) > v_cc + 1e-9))
    if not beyond:
        return ""
    return (
        f"{beyond} of {x.size} points lie beyond |V| = {v_cc:.0f} V, the deepest "
        "bias at which the batched engine is cross-checked against devsim. devsim "
        "itself may not converge there (default device: about −64 V), so these "
        "values are not validated."
    )


def truncation_text(result) -> str:
    md = result.metadata
    if not md.get("truncated") or len(result.x) == 0:
        return ""
    return (
        f"Sweep stopped at {np.min(result.x):.1f} V (requested down to "
        f"{md.get('requested_v_stop', float('nan')):.1f} V): devsim did not "
        "converge beyond this point. This is an artefact of the device model "
        "(ohmic cathode directly on the n⁻ epi, no n⁺ substrate) and of solver "
        "robustness beyond full depletion (outside the validated domain), not a "
        "physical device limit. Results below cover only "
        "the converged portion."
    )


def show_result_notes(result) -> None:
    """Solver badge, truncation warning and full-depletion caption for one result."""
    trunc = truncation_text(result)
    if trunc:
        st.warning(trunc)
    badge = solver_badge_text(result.metadata)
    if badge:
        st.caption(badge)
    rt = reach_through_text(result)
    if rt:
        st.caption(rt)
    cc = cross_check_text(result)
    if cc:
        st.warning(cc)

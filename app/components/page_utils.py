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

from app.components.results import provenance_text

SOLVER_LABELS = {
    "batched": "batched 1D engine (numba)",
    "engine2d": "2D engine (numba)",
    "batched+traps": "batched 1D engine with charged deep traps",
    "devsim": "devsim (reference finite-volume solver)",
}

STALE_RESULT_TEXT = (
    "These results are for a previous configuration. The sidebar settings have "
    "changed since this run; press **Run simulation** to update."
)
STALE_INPUTS_TEXT = (
    "These results were computed with different page inputs ({changed}); press "
    "**Run simulation** to update."
)
EMPTY_CONFIG_TEXT = "Configure a device in the sidebar to begin."
INVALID_CONFIG_TEXT = (
    "The device configuration in the sidebar is invalid, so nothing can be run: "
    "{error}. Fix the value in the sidebar."
)
CROSS_CHECK_SCOPE = (
    "cross-check = the batched engine and devsim agree (relative difference "
    "≤ 3e-6, typically far smaller) at the checked biases down to −40 V "
    "(−100 V with the n+ layer); neither is validated against measurement there"
)


def about(text: str) -> None:
    """One collapsed 'About this page' box: what the page computes, the model, its
    validity limits (Markdown)."""
    with st.expander("About this page: model and validity", expanded=False):
        st.markdown(text)


def _cathode_text(result) -> str:
    """The cathode model the saturation / truncation depends on, per the result's config."""
    cfg = getattr(result, "config", None)
    if cfg is not None and getattr(cfg, "nplus_thickness_um", None) is None:
        return "ohmic contact directly on the n⁻ epi, no n⁺ layer (legacy model)"
    return "the 1 µm n⁺ layer under the epi, cathode on its far end"


def require_config():
    """The current DeviceConfig, or stop the page with a clear message.

    Never falls back to a stale valid config: when the sidebar holds an invalid
    value, ``device_config`` is absent and ``device_config_error`` explains why.
    """
    cfg = st.session_state.get("device_config")
    if cfg is not None:
        return cfg
    err = st.session_state.get("device_config_error")
    if err:
        st.error(INVALID_CONFIG_TEXT.format(error=err))
    else:
        st.info(EMPTY_CONFIG_TEXT)
    st.stop()


def show_error(exc: BaseException, hint: str = "") -> None:
    """Friendly st.error plus an expander with the technical details."""
    name = type(exc).__name__
    if name == "OutsideModelValidity":  # message written for the user
        head = str(exc)
    elif isinstance(exc, ValueError):
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


def show_stale_warning(result_config, cfg, used_inputs=None, current_inputs=None) -> None:
    """Warn when the cached result was computed for a different sidebar config,
    or (when given) with different page inputs than the ones now on the page."""
    if result_config is not None and cfg is not None and result_config != cfg:
        st.warning(STALE_RESULT_TEXT)
    if used_inputs is not None and current_inputs is not None:
        changed = [
            k for k in current_inputs if used_inputs.get(k) != current_inputs[k]
        ]
        if changed:
            st.warning(STALE_INPUTS_TEXT.format(changed=", ".join(changed)))


def solver_badge_text(metadata: dict) -> str:
    """'Solver: ... · inside/outside the engine/devsim cross-check domain'."""
    solver = metadata.get("solver")
    if solver is None:
        return ""
    # mixed results name several solvers ("batched + devsim")
    label = " + ".join(SOLVER_LABELS.get(s, s) for s in str(solver).split(" + "))
    parts = [f"Solver: {label}"]
    if metadata.get("fallback_reason"):
        parts[0] += f" (fallback: {metadata['fallback_reason']})"
    text = provenance_text(metadata)
    if text:
        parts.append(text)
    return " · ".join(parts)


_SATURATES = {
    "cv": "the capacitance saturates",
    "cce": "the collected charge saturates",
}


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
            "Full-depletion voltage not computed (device outside the "
            "engine/devsim cross-check domain)."
        )
    what = _SATURATES.get(result.sim_type, "the result saturates")
    method = md.get("reach_through_method")
    return (
        f"Epi fully depleted from |V| ≈ {v_rt:.2f} V (dotted line"
        + (f"; {method}" if method else "")
        + f"). Beyond it {what} and the result depends on the cathode model "
        f"({_cathode_text(result)})."
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
        f"{beyond} of {x.size} {'point lies' if beyond == 1 else 'points lie'} "
        f"beyond |V| = {v_cc:.0f} V, the deepest "
        "bias at which the batched engine is cross-checked against devsim. devsim "
        "itself may not converge there (legacy model without n+: about −64 V), so "
        "these values are not validated."
    )


def truncation_text(result) -> str:
    md = result.metadata
    if not md.get("truncated") or len(result.x) == 0:
        return ""
    return (
        f"Sweep stopped at {np.min(result.x):.1f} V (requested down to "
        f"{md.get('requested_v_stop', float('nan')):.1f} V): devsim did not "
        "converge beyond this point. This is an artefact of the device model "
        f"({_cathode_text(result)}) and of solver "
        "robustness beyond full depletion (outside the cross-check domain), not a "
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

"""UI-01 coverage: empty-state guard unit test + AppTest nav/sidebar smoke.

test_empty_state_guard uses AppTest.from_function on a small wrapper that
imports and calls app.workflows.cv.render() — AppTest.from_function requires
the supplied callable's body to be self-contained (its own imports), so the
actual page render() is imported *inside* the wrapper body rather than at
module scope. This still exercises the real guard logic defined in
app/workflows/cv.py, not a copy of it.

test_nav_sidebar_smoke boots the full app/main.py entry script headlessly.
"""

from pathlib import Path

from streamlit.testing.v1 import AppTest


def test_first_edit_survives_rerun():
    """UI-07 first-rerun persistence regression (VERIFICATION.md gap).

    Boots the REAL app/main.py, sets the Epi-thickness sidebar widget to a
    non-default 7.5 on script-run #1, reruns once (the first post-boot rerun),
    and asserts device_config.epi_thickness_um is still 7.5. Against the
    current unmodified app this FAILS: the value silently reverts to the
    hardcoded 10.0 default because the app/pages/ magic-dir collision with
    st.navigation destabilizes implicit widget-ID computation across the
    run1->run2 transition (38-VERIFICATION.md, Truth 3).

    Accessor note: at.sidebar.number_input[0] addresses "Epi thickness (µm)"
    in the default graded/1D render (empirically confirmed on Streamlit
    1.55.0 — it is the first number_input rendered; set_value(7.5)+run()
    round-trips).
    """
    at = AppTest.from_file("app/main.py")
    at.run()  # script-run #1 (initial boot)

    epi_widget = at.sidebar.number_input[0]
    assert epi_widget.label == "Epi thickness (µm)", (
        f"accessor drift: sidebar number_input[0] is {epi_widget.label!r}, "
        "expected 'Epi thickness (µm)'"
    )
    epi_widget.set_value(7.5)

    at.run()  # script-run #2 (first rerun)

    assert at.session_state["device_config"].epi_thickness_um == 7.5, (
        "first post-boot sidebar edit was discarded on the first rerun — "
        "reverting to the 10.0 default is the VERIFICATION.md UI-07 gap "
        f"(got {at.session_state['device_config'].epi_thickness_um})"
    )


def test_no_magic_pages_directory():
    """Structural root-cause guard against a sibling app/pages directory.

    The directory rename (app/pages -> app/workflows) is the empirically
    load-bearing fix against the real app: adding explicit key= alone, without
    the rename, does not close the defect (38-VERIFICATION.md). This guard
    checks the root cause DIRECTLY — a future edit could reintroduce a sibling
    app/pages directory without any behavioral test catching it (a behavioral
    test could pass for an unrelated reason), so we assert the pages-named
    directory never returns and that app/workflows exists in its place.
    """
    repo_root = Path(__file__).resolve().parents[1]
    assert not (repo_root / "app" / "pages").exists(), (
        "app/pages exists — a sibling directory literally named 'pages' next "
        "to app/main.py triggers Streamlit legacy magic-multipage auto-detection "
        "colliding with st.navigation (the VERIFICATION.md root cause). It must "
        "be renamed to app/workflows and never reintroduced."
    )
    workflows = repo_root / "app" / "workflows"
    assert workflows.is_dir(), (
        "app/workflows does not exist as a directory — the page modules must "
        "live under app/workflows (renamed from app/pages)."
    )


def test_empty_state_guard():
    def _run_cv_page():
        from app.workflows.cv import render

        render()

    at = AppTest.from_function(_run_cv_page)
    # session_state is empty by default (do not pre-seed device_config).
    at.run()

    assert at.exception == [], f"page crashed on empty session_state: {at.exception}"
    info_texts = [el.value for el in at.info]
    assert "Configure a device in the sidebar to begin." in info_texts


def test_nav_sidebar_smoke():
    at = AppTest.from_file("app/main.py")
    at.run()

    # AppTest's public API does not expose a nav-page-count accessor (no
    # element type for st.navigation's sidebar entries), so this asserts
    # the documented minimum smoke: the entry script (which registers all
    # 8 pages via st.navigation before pg.run()) boots with no exception,
    # and the sidebar wrote device_config before running the default page.
    assert at.exception == [], f"app/main.py raised on boot: {at.exception}"
    assert "device_config" in at.session_state
    assert at.session_state["device_config"] is not None

"""UI-01 coverage: empty-state guard unit test + AppTest nav/sidebar smoke.

test_empty_state_guard uses AppTest.from_function on a small wrapper that
imports and calls app.pages.cv.render() — AppTest.from_function requires the
supplied callable's body to be self-contained (its own imports), so the
actual page render() is imported *inside* the wrapper body rather than at
module scope. This still exercises the real guard logic defined in
app/pages/cv.py, not a copy of it.

test_nav_sidebar_smoke boots the full app/main.py entry script headlessly.
"""

from streamlit.testing.v1 import AppTest


def test_empty_state_guard():
    def _run_cv_page():
        from app.pages.cv import render

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

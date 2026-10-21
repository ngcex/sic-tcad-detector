"""Smoke test of an INSTALLED etna (wheel), run outside the source checkout.

    uv build -o dist
    uv venv /tmp/v -p 3.13 && VIRTUAL_ENV=/tmp/v uv pip install dist/etna-*.whl gdstk
    cd /tmp && /tmp/v/bin/python -I /path/to/checkout/scripts/wheel_smoke.py [--app]

Fails (non-zero exit) when etna is imported from anywhere but site-packages, when the
packaged data (etna/_data) is not the one in use, or when any step raises. ``--app``
also copies ``app/`` (not part of the wheel) to a temp dir, renders the entry script and
the pages that read packaged data with streamlit's AppTest, and checks the server's
``/_stcore/health`` endpoint. Used by the ``wheel-smoke`` CI job.
"""

from __future__ import annotations

import argparse
import shutil
import socket
import subprocess
import sys
import sysconfig
import tempfile
import time
import urllib.request
import warnings
from pathlib import Path

REPO = Path(__file__).resolve().parents[1]  # only to copy app/ for --app


def step(name, fn):
    t0 = time.perf_counter()
    out = fn()
    print(f"[ok] {name} ({time.perf_counter() - t0:.1f} s){': ' + out if out else ''}")


def check_install():
    import etna
    from etna import _paths

    site = Path(sysconfig.get_paths()["purelib"]).resolve()
    where = Path(etna.__file__).resolve()
    assert site in where.parents, f"etna imported from {where}, not {site}"
    root = _paths.data_root().resolve()
    assert root == where.parent / "_data", f"data root {root} is not the packaged one"
    return f"etna {etna.__version__} at {where.parent}, data {root}"


def stopping_power():
    from etna.core import stopping_power as sp

    t = sp.load_star("proton", "water")
    assert len(t["energy_MeV"]) > 50
    return f"PSTAR water {len(t['energy_MeV'])} rows"


def kappa():
    from etna.core.microdosimetry import compute_kappa_table

    k = compute_kappa_table(source="bragg")["kappa"]
    assert 0.3 < float(k.min()) < float(k.max()) < 0.5, (k.min(), k.max())
    with warnings.catch_warnings():
        warnings.simplefilter("ignore")
        legacy = compute_kappa_table(source="legacy_unsafe")["kappa"]
    assert len(legacy) > 0
    return f"bragg kappa {float(k.min()):.3f}-{float(k.max()):.3f}; legacy table read"


def run_cv():
    import etna

    r = etna.run_cv(etna.DeviceConfig())
    assert len(r.x) > 2 and not r.metadata.get("error"), r.metadata.get("error")
    return f"{len(r.x)} points"


def run_cce():
    import etna

    r = etna.run_cce(etna.DeviceConfig(), n_points=5)
    assert len(r.y) == 5 and 0.0 <= float(max(r.y)) <= 1.0 + 1e-9, r.y
    return f"CCE {float(min(r.y)):.3f}-{float(max(r.y)):.3f}"


def run_all():
    import etna

    rep = etna.run_all(etna.DeviceConfig())
    assert rep.ok, rep.rows()
    return f"{len(rep.sections)} sections ok"


def viz_validation():
    import etna.viz.validation as v

    r = v.w_creep()
    assert r.metadata, "w_creep metadata empty"
    return "w_creep ok"


def layout():
    import etna.export as ex

    p = ex.load_layout_params()
    assert Path(ex.LAYOUT_PARAMS_PATH).is_file()
    with warnings.catch_warnings():
        warnings.simplefilter("ignore", ex.PlaceholderDimensionWarning)
        svg = ex.layout_svg("dosimeter", p)
        assert svg.lstrip().startswith("<")
        ok, why = ex.gds_available()
        if ok:
            gds = ex.layout_gds("dosimeter", p)
            assert len(gds) > 100
            extra = f", GDS {len(gds)} B"
        else:
            extra = f", GDS skipped ({why})"
        cfg = __import__("etna").DeviceConfig()
        spec = ex.build_spec(cfg, layout_params=None)
        pdf = ex.spec_to_pdf(spec, config=cfg)  # re-reads the params file from the spec
        assert pdf[:4] == b"%PDF"
    return (
        f"params {p.display_path}, SVG {len(svg)} chars{extra}, spec PDF {len(pdf)} B"
    )


def app_checks():
    from streamlit.testing.v1 import AppTest

    import etna

    tmp = Path(tempfile.mkdtemp(prefix="etna-app-"))
    shutil.copytree(REPO / "app", tmp / "app")
    sys.path.insert(0, str(tmp))
    import app.workflows.figures as fig
    import app.workflows.measurements as meas

    for name, path in fig.LET_TABLES.items():
        etna.viz.load_let_table(path)
    templates = sorted(p.name for p in meas.TEMPLATES.glob("*.csv"))
    assert len(templates) == 3, templates

    at = AppTest.from_file(str(tmp / "app" / "main.py"), default_timeout=240)
    at.run()
    assert at.exception == [], at.exception
    mod = sys.modules["etna"].__file__
    assert "site-packages" in mod, mod

    def _measurements():
        from app.workflows.measurements import render

        render()

    at = AppTest.from_function(_measurements, default_timeout=240)
    at.session_state["device_config"] = etna.DeviceConfig()
    at.run()
    assert at.exception == [], at.exception
    n_dl = len(at.get("download_button"))
    assert n_dl >= 3, f"{n_dl} template download buttons"

    with socket.socket() as s:
        s.bind(("127.0.0.1", 0))
        port = s.getsockname()[1]
    cmd = [
        sys.executable,
        "-m",
        "streamlit",
        "run",
        str(tmp / "app" / "main.py"),
        "--server.headless",
        "true",
        "--server.address",
        "127.0.0.1",
        "--server.port",
        str(port),
        "--browser.gatherUsageStats",
        "false",
    ]
    proc = subprocess.Popen(
        cmd, cwd=tmp, stdout=subprocess.PIPE, stderr=subprocess.STDOUT
    )
    try:
        url = f"http://127.0.0.1:{port}/_stcore/health"
        for _ in range(120):
            try:
                body = urllib.request.urlopen(url, timeout=2).read().decode()
                break
            except OSError:
                time.sleep(0.5)
        else:
            raise AssertionError("streamlit server never answered /_stcore/health")
        assert body.strip() == "ok", body
    finally:
        proc.terminate()
        proc.wait(timeout=30)
    shutil.rmtree(tmp, ignore_errors=True)
    return (
        f"LET tables {list(fig.LET_TABLES)}, templates {templates}, "
        f"{n_dl} download buttons, health ok"
    )


def main():
    ap = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    ap.add_argument("--app", action="store_true", help="also smoke the Streamlit app")
    ap.add_argument("--quick", action="store_true", help="skip run_all")
    args = ap.parse_args()
    if (Path.cwd() / "etna").is_dir():
        sys.exit("run from outside the source checkout (cwd has an etna/ directory)")
    step("installed package + packaged data", check_install)
    step("NIST stopping power", stopping_power)
    step("SiC->water kappa (bragg + legacy tables)", kappa)
    step("run_cv(DeviceConfig())", run_cv)
    step("run_cce(DeviceConfig())", run_cce)
    if not args.quick:
        step("run_all(DeviceConfig())", run_all)
    step("viz.validation.w_creep", viz_validation)
    step("layout params, SVG/GDS, spec PDF", layout)
    if args.app:
        step("Streamlit app", app_checks)
    print("wheel smoke: all passed")


if __name__ == "__main__":
    main()

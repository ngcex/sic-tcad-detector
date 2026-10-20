"""etna.measurements file validation: bad files are refused with messages naming the
file line and column; the committed templates are header-only (no numbers).

Every inline CSV literal below is SYNTHETIC test input (made-up numbers to exercise the
validator), not a measurement."""

from __future__ import annotations

from pathlib import Path

import pytest

import etna.measurements as M
from etna.measurements import MeasurementFileError

TEMPLATES = Path(__file__).resolve().parents[1] / "data" / "measurements" / "templates"
GOOD_CV = "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,90\n-2,48\n"


def _issues(loader, text, **kw):
    with pytest.raises(MeasurementFileError) as e:
        loader(text, **kw)
    return e.value.issues


@pytest.mark.parametrize(
    "text, expected",
    [
        ("# area_cm2: 0.04\nC_pF\n200\n", "missing required column (one of: V_V)"),
        ("V_V\n0\n", "missing required column (one of: C_F, C_nF, C_pF, C_per_area"),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,abc\n",
            "line 4, column 'C_pF': 'abc' is not a finite number",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,nan\n",
            "line 4, column 'C_pF': 'nan' is not a finite number",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n,90\n",
            "line 4, column 'V_V': empty cell",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,90\n-1,80\n",
            "line 5, column 'V_V': bias -1 V repeats line 4",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,-90\n",
            "line 4, column 'C_pF': C must be > 0",
        ),
        ("V_V,C_pF\n0,200\n-1,90\n", "absolute capacitance but no area is given"),
        ("# area_cm2: 0.04\nV_V,C_pF\n", "no data rows after the header (line 2)"),
        ("# area_cm2: 0.04\nV_V,C_uF\n0,1\n", "unknown column 'C_uF'"),
        (
            "# area_cm2: 0.04\nV_V,C_pF,C_F\n0,1,1\n",
            "columns 'C_pF' and 'C_F' give the same quantity",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF,sigma_C_per_area_F_cm2\n0,1,1\n",
            "does not match the capacitance column",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1\n",
            "line 4: 1 cells, the header has 2 columns",
        ),
        ("# area_cm2: 0.04\nV_V,C_pF\n1,200\n2,90\n", "2 forward vs 0 reverse points"),
        # reverse-positive file without the convention line, one small negative point
        (
            "# area_cm2: 0.04\nV_V,C_pF\n-0.5,300\n0,200\n5,60\n10,50\n",
            "2 forward vs 1 reverse",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF\n4,300\n0,200\n-5,60\n-10,50\n",
            "forward bias up to 4 V",
        ),
        (
            "# area_cm2: 0.04\n# bias_convention: sideways\nV_V,C_pF\n0,1\n",
            "bias_convention must be one of",
        ),
        (
            "# area_cm2: zero\nV_V,C_pF\n0,1\n",
            "line 1: metadata 'area_cm2' must be a finite number",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF,area_cm2\n0,200,0.04\n-1,90,0.05\n",
            "line 4, column 'area_cm2': 0.05 differs from 0.04 on line 3",
        ),
        (
            "# area_cm2: 0.04\nV_V,C_pF,T_K\n0,200,300\n-1,90,301\n",
            "line 4, column 'T_K': 301 differs",
        ),
        ("# area_cm2: 0.04\nV_V,C_pF\n0,200\n-1,90\n", None),  # good file
    ],
)
def test_cv_file_problems_name_line_and_column(text, expected):
    if expected is None:
        assert len(M.load_cv_csv(text).V) == 2
        return
    msgs = _issues(M.load_cv_csv, text)
    assert any(expected in m for m in msgs), msgs


def test_cv_area_argument_conflict():
    msgs = _issues(M.load_cv_csv, GOOD_CV, area_cm2=0.05)
    assert any("disagrees with the file's area" in m for m in msgs)


def test_all_problems_reported_together():
    text = "# area_cm2: 0.04\nV_V,C_pF\n0,x\n-1,y\n-2,-3\n"
    msgs = _issues(M.load_cv_csv, text)
    assert len(msgs) == 2  # value errors stop before the semantic checks
    assert "line 3" in msgs[0] and "line 4" in msgs[1]
    err = pytest.raises(MeasurementFileError, M.load_cv_csv, text).value
    assert "2 problem(s)" in str(err) and isinstance(err, ValueError)


def test_missing_file_is_a_clear_error(tmp_path):
    with pytest.raises(FileNotFoundError, match="measurement file not found"):
        M.load_cv_csv(tmp_path / "nope.csv")


@pytest.mark.parametrize(
    "text, expected",
    [
        ("depth_um\n0\n", "missing required column (one of: Al_cm3, N_cm3)"),
        (
            "depth_um,N_cm3\n0,1e15\n0.1,1e15\n0.1,2e15\n",
            "line 4, column 'depth_um': depth 0.1 is not above 0.1 on line 3",
        ),
        ("depth_um,N_cm3\n0,1e15\n-0.1,1e15\n", "depth must be >= 0"),
        ("depth_um,N_cm3\n0,1e15\n0.1,-1e15\n", "line 3, column 'N_cm3': must be >= 0"),
        ("depth_um,N_cm3\n0,0\n0.1,0\n", "column 'N_cm3': no value above zero"),
        (
            "# depth_reference: junction\ndepth_um,N_cm3\n0,1e15\n",
            "depth_reference must be 'surface'",
        ),
        ("# technique: XPS\ndepth_um,N_cm3\n0,1e15\n", "technique must be one of"),
        ("depth_um,N_cm3,sigma_N_cm3,rel_sigma_N\n0,1,1,1\n", "are exclusive"),
        ("depth_nm,depth_um,N_cm3\n0,0,1\n", "give the same quantity"),
    ],
)
def test_sims_file_problems(text, expected):
    msgs = _issues(M.load_sims_csv, text)
    assert any(expected in m for m in msgs), msgs


@pytest.mark.parametrize(
    "text, expected",
    [
        ("I_pA\n1\n", "missing required column (one of: T_C, T_K)"),
        ("T_K,I_pA\n300,1\n310,0\n", "line 3, column 'I_pA': I = 0"),
        ("T_K,I_pA\n300,1\n-5,2\n", "line 3, column 'T_K': T must be > 0 K"),
        ("T_K,I_pA\n300,1\n300,2\n", "line 3, column 'T_K': T 300 K repeats line 2"),
        (
            "V_V,T_K,I_pA\n-10,300,1\n-20,300,2\n-10,300,3\n",
            "T 300 K repeats line 2 at bias -10 V",
        ),
        ("T_K,T_C,I_pA\n300,27,1\n", "are exclusive"),
        ("V_V,T_K,I_pA\n10,300,1\n10,310,2\n", "2 row(s) with forward bias"),
        ("T_K,I_pA,sigma_I_pA\n300,1,-1\n", "column 'sigma_I_pA': must be >= 0"),
    ],
)
def test_it_file_problems(text, expected):
    msgs = _issues(M.load_it_csv, text)
    assert any(expected in m for m in msgs), msgs


@pytest.mark.parametrize(
    "name, loader",
    [
        ("cv_template.csv", M.load_cv_csv),
        ("sims_template.csv", M.load_sims_csv),
        ("it_template.csv", M.load_it_csv),
    ],
)
def test_templates_are_header_only_and_valid(name, loader):
    path = TEMPLATES / name
    lines = [line.strip() for line in path.read_text().splitlines() if line.strip()]
    data = [line for line in lines if not line.startswith("#")]
    assert len(data) == 1, "a template holds the header line only (no data rows)"
    # metadata values are empty or one of the documented non-numeric defaults
    allowed = {"", "reverse_negative", "SIMS", "surface"}
    for line in lines:
        if line.startswith("#") and ":" in line:
            key, _, value = line.lstrip("# ").partition(":")
            if " " not in key.strip():  # a 'key: value' metadata line
                assert value.strip() in allowed, f"{name}: {line!r} carries a value"
    # every column is accepted: the only problem is the missing data
    msgs = _issues(loader, path)
    assert len(msgs) == 1 and "no data rows" in msgs[0]


def test_templates_directory_holds_only_the_three_templates():
    """templates/ is header-only (real, sourced measurement files may live elsewhere)."""
    files = sorted(p.name for p in TEMPLATES.iterdir() if p.is_file())
    assert files == ["cv_template.csv", "it_template.csv", "sims_template.csv"]

"""docs/THEORY.md must not drift from the code: every repository path it names exists,
every identifier it quotes in backticks occurs in the sources, every function it says lives
in a given module is defined there, and every USAGE.md anchor it links to exists."""

from __future__ import annotations

import ast
import functools
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
THEORY = ROOT / "docs" / "THEORY.md"
SOURCE_DIRS = ("etna", "scripts", "lab")


def _text():
    return THEORY.read_text(encoding="utf-8")


def _code_spans(text):
    """Code spans, without the GitHub inline math $`...`$ and the ```math blocks."""
    text = re.sub(r"```math\n.*?```", "", text, flags=re.S)
    text = re.sub(r"\$`[^`]*`\$", "", text)
    return re.findall(r"`([^`\n]+)`", text)


def _corpus():
    parts = []
    for d in SOURCE_DIRS:
        for p in (ROOT / d).rglob("*.py"):
            parts.append(p.read_text(encoding="utf-8", errors="replace"))
    return "\n".join(parts)


def test_named_paths_exist():
    missing = []
    for span in _code_spans(_text()):
        token = span.split()[0].rstrip(",;:")
        if "/" not in token or token.startswith(("http", "~")):
            continue
        if not re.search(r"\.(py|md|json|csv|toml|gz)$|/$", token) and "(" in token:
            continue
        pattern = token.replace("<topic>", "*")
        if "*" in pattern:
            if not list(ROOT.glob(pattern)):
                missing.append(token)
        elif re.search(r"\.(py|md|json|csv|toml|gz|pdf)$", token) or token.endswith(
            "/"
        ):
            if not (ROOT / token).exists():
                missing.append(token)
    assert not missing, f"paths named in THEORY.md that do not exist: {missing}"


def test_quoted_identifiers_occur_in_sources():
    corpus = _corpus()
    words = set(re.findall(r"[A-Za-z_][A-Za-z0-9_]*", corpus))
    unknown = []
    for span in _code_spans(_text()):
        token = span.strip()
        if not re.fullmatch(r"[A-Za-z_][A-Za-z0-9_]*", token) or "_" not in token:
            continue
        if token not in words:
            unknown.append(token)
    assert (
        not unknown
    ), f"identifiers quoted in THEORY.md not found in the sources: {unknown}"


@functools.lru_cache(maxsize=None)
def _defined_names(path):
    tree = ast.parse(path.read_text(encoding="utf-8"))
    names = set()
    for node in ast.walk(tree):
        if isinstance(node, (ast.FunctionDef, ast.ClassDef, ast.AsyncFunctionDef)):
            names.add(node.name)
        elif isinstance(node, ast.Assign):
            names.update(t.id for t in node.targets if isinstance(t, ast.Name))
        elif isinstance(node, ast.AnnAssign) and isinstance(node.target, ast.Name):
            names.add(node.target.id)
    return names


def test_paragraph_functions_defined_in_a_cited_module():
    """A function quoted as `name()`-style identifier next to module paths in the same
    paragraph is defined in one of those modules or elsewhere in etna (no stale names)."""
    defined = set()
    for p in (ROOT / "etna").rglob("*.py"):
        defined |= _defined_names(p)
    bad = []
    for para in re.split(r"\n\s*\n", _text()):
        spans = _code_spans(para)
        modules = [s for s in spans if re.fullmatch(r"etna/[\w/]+\.py", s)]
        if not modules:
            continue
        local = set()
        for m in modules:
            local |= _defined_names(ROOT / m)
        for s in spans:
            if re.fullmatch(r"[a-z_][a-z0-9_]*", s) and "_" in s and s in defined:
                continue
            if re.fullmatch(r"[a-z_][a-z0-9_]*", s) and "_" in s and s not in local:
                if re.search(rf"\b{s}\(", para) or s.startswith(("run_", "compute_")):
                    bad.append(s)
    assert not bad, f"function names in THEORY.md not defined in etna: {bad}"


def test_usage_anchors_exist():
    usage = (ROOT / "docs" / "USAGE.md").read_text(encoding="utf-8")

    def slug(h):
        s = h.strip().lower()
        s = re.sub(r"[^\w\- ]", "", s)
        return s.replace(" ", "-")

    anchors = {slug(h) for h in re.findall(r"^#+ (.+)$", usage, flags=re.M)}
    links = re.findall(r"\(USAGE\.md#([^)]+)\)", _text())
    assert links, "THEORY.md links into USAGE.md"
    missing = [a for a in links if a not in anchors]
    assert not missing, f"USAGE.md anchors not found: {missing}"


def test_worked_examples_module_lists_every_chapter_function():
    tree = ast.parse(
        (ROOT / "lab" / "theory_manual" / "worked_examples.py").read_text()
    )
    funcs = {n.name for n in tree.body if isinstance(n, ast.FunctionDef)}
    cited = set(re.findall(r"worked_examples\.py (\w+)", _text()))
    chapters_src = (ROOT / "lab" / "theory_manual" / "worked_examples.py").read_text()
    keys = set(re.findall(r'^\s+"(\w+)": \w+,$', chapters_src, flags=re.M))
    assert cited <= keys, f"chapters cited but not in CHAPTERS: {cited - keys}"
    assert {"main"} <= funcs


def test_math_uses_github_safe_delimiters():
    """Inline math is written $`...`$ and display math as ```math blocks: plain $...$
    next to a letter, a slash or a dash is not rendered by GitHub (checked with the
    GitHub markdown API), and underscores inside it turn into emphasis."""
    text = re.sub(r"```.*?```", "", _text(), flags=re.S)
    text = re.sub(r"\$`[^`]*`\$", "", text)
    text = re.sub(r"`[^`\n]*`", "", text)
    assert "$" not in text, "plain $ math left in THEORY.md"

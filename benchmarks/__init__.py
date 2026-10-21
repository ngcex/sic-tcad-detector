"""ETNA regression-benchmark registry (issue #3, P2.2).

``registry.toml`` lists the golden quantities (id, quantity, method, reference with
uncertainty and source, tolerance, status class); ``evaluators.py`` computes them;
``runner.py`` compares and writes JSON + Markdown. Run from the repository root:

    uv run python -m benchmarks            # fast subset -> benchmarks/results/fast.{json,md}
    uv run python -m benchmarks --slow     # everything evaluable -> results/all.*
    uv run pytest -m benchmark -q          # the same fast subset as tests

Not shipped in the wheel (a development tool, like tests/). The V&V report is
docs/VV_REPORT.md.
"""

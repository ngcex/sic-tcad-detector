"""LIB-05 vertical-slice validation script.

Proves the DeviceConfig -> build_device -> core -> SimResult contract works
end-to-end through the public `etna` API: build a default device
configuration, run a C-V sweep over the facade's default range (0 to -50 V),
print the bias and capacitance arrays, and check two things, raising if either
fails:

- the sweep was not truncated (it reached the requested end bias), and
- the capacitance decreases monotonically as reverse bias increases (more
  negative bias -> lower C).

Run: python examples/cv_example.py
(tests/test_examples.py runs it as a smoke test.)
"""

from __future__ import annotations

import numpy as np

from etna import DeviceConfig, run_cv


def main() -> None:
    cfg = DeviceConfig()  # all defaults -> 1D graded device

    result = run_cv(cfg)  # facade default: 0 .. -50 V

    print("Bias (V):")
    print(result.x)
    print()
    print("Capacitance (F):")
    print(result.y)
    print()

    truncated = bool(result.metadata.get("truncated", False))
    requested = result.metadata.get("requested_v_stop")
    print(
        f"truncated={truncated}  requested_v_stop={requested}  "
        f"last bias reached={result.x[-1]:.1f} V  solver={result.metadata.get('solver')}"
    )
    if truncated:
        raise RuntimeError(
            f"C-V sweep truncated at {result.x[-1]:.1f} V "
            f"(requested down to {requested} V)"
        )

    # Order from 0 V towards the deepest reverse bias, then require dC/d|V| < 0.
    order = np.argsort(-np.asarray(result.x))
    C = np.asarray(result.y)[order]
    if not np.all(np.diff(C) < 0):
        raise RuntimeError(
            "capacitance does not decrease monotonically with reverse bias"
        )
    print(
        f"C at {result.x[order[0]]:.1f}V = {C[0]:.4e} F, "
        f"C at {result.x[order[-1]]:.1f}V = {C[-1]:.4e} F "
        "(decreasing monotonically with reverse bias)"
    )


if __name__ == "__main__":
    main()

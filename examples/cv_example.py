"""LIB-05 vertical-slice validation script.

Proves the DeviceConfig -> build_device -> core -> SimResult contract works
end-to-end through the public `etna` API: build a default device
configuration, run a C-V sweep, and print the resulting bias and
capacitance arrays. Capacitance must decrease monotonically as reverse bias
increases (more negative bias -> lower C).

Run: python examples/cv_example.py
"""

from __future__ import annotations

from etna import DeviceConfig, run_cv


def main() -> None:
    cfg = DeviceConfig()  # all defaults -> 1D graded device

    result = run_cv(cfg, v_start=0, v_stop=-200, n_points=20)

    print("Bias (V):")
    print(result.x)
    print()
    print("Capacitance (F):")
    print(result.y)
    print()

    c_at_0v = result.y[0]
    c_at_most_reverse = result.y[-1]
    print(
        f"C at {result.x[0]:.1f}V = {c_at_0v:.4e} F, "
        f"C at {result.x[-1]:.1f}V = {c_at_most_reverse:.4e} F "
        f"({'decreasing' if c_at_most_reverse < c_at_0v else 'NOT decreasing'} "
        "with reverse bias)"
    )


if __name__ == "__main__":
    main()

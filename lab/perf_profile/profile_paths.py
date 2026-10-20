"""Profile user-facing ETNA paths on the default config, one subprocess per path.

usage:
  uv run python lab/perf_profile/profile_paths.py            # all paths
  uv run python lab/perf_profile/profile_paths.py run_cv cce  # substring filter
  uv run python lab/perf_profile/profile_paths.py --top 25 run_field_1d

Per path: warm-up call (numba cache, imports), then one cProfile'd call.
Reports wall, cpu (process time, all threads), and tottime shares by bucket:
devsim (devsim.* builtins), numba (dispatch/compile + jitted kernels, seen as the
dispatcher call), numpy, scipy, pandas, matplotlib, io, python (everything else, i.e.
etna python loops). Results append to lab/perf_profile/results/paths.jsonl.

NOTE: the machine is shared with other agents (load ~30-45 on 10 cores at the time of
writing), so absolute wall times are inflated; cpu and the shares are more robust.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
import time
from pathlib import Path

HERE = Path(__file__).resolve().parent
OUT = HERE / "results"

# name -> (setup code, call expression); executed in a fresh interpreter.
PATHS: dict[str, tuple[str, str]] = {
    "run_cv": ("from etna import DeviceConfig, run_cv; c=DeviceConfig()", "run_cv(c)"),
    "run_cce": (
        "from etna import DeviceConfig, run_cce; c=DeviceConfig()",
        "run_cce(c)",
    ),
    "run_field_1d": (
        "from etna import DeviceConfig, run_field; c=DeviceConfig()",
        "run_field(c, -30.0)",
    ),
    "run_field_2d": (
        "from etna import DeviceConfig, run_field; c=DeviceConfig(half_width_um=50.0)",
        "run_field(c, -30.0)",
    ),
    "run_temperature_sweep_dd": (
        "from etna import DeviceConfig, run_temperature_sweep; c=DeviceConfig()",
        "run_temperature_sweep(c, method='dd')",
    ),
    "run_dark_current": (
        "from etna import DeviceConfig, run_dark_current; c=DeviceConfig()",
        "run_dark_current(c)",
    ),
    "run_dark_current_devsim": (
        "from etna import DeviceConfig, run_dark_current; c=DeviceConfig()",
        "run_dark_current(c, 0.0, -60.0, 7, solver='devsim')",
    ),
    "run_flash": (
        "from etna import DeviceConfig, run_flash_recombination; c=DeviceConfig()",
        "run_flash_recombination(c)",
    ),
    "run_transient": (
        "from etna import DeviceConfig, run_transient; c=DeviceConfig()",
        "run_transient(c)",
    ),
    "run_radiation_damage": (
        "from etna import DeviceConfig, run_radiation_damage; c=DeviceConfig()",
        "run_radiation_damage(c)",
    ),
    "run_microdosimetry": (
        "from etna import DeviceConfig, run_microdosimetry; c=DeviceConfig()",
        "run_microdosimetry(c, 'data/synthetic_mc_events.csv')",
    ),
    "run_cce_tilted_alpha": (
        "from etna import microdosimeter_config, run_cce_tilted_alpha; c=microdosimeter_config()",
        "run_cce_tilted_alpha(c)",
    ),
    "run_all": (
        "from etna import DeviceConfig, run_all; c=DeviceConfig()",
        "run_all(c)",
    ),
    "run_all_sweep": (
        "from etna import DeviceConfig, run_all; c=DeviceConfig()",
        "run_all(c, sections=('sweep',))",
    ),
    "dark_current_2d": (
        "from etna.core.optimization import get_dark_current_2d",
        "get_dark_current_2d(half_width_um=50.0, V_bias=-40.0)",
    ),
    "cce_uncertainty_envelope": (
        "import numpy as np\nfrom etna.core.radiation_damage import cce_uncertainty_envelope",
        "cce_uncertainty_envelope(np.geomspace(1e11, 1e13, 6), V_bias=-40)",
    ),
    "radiation_hardness_sweep": (
        "from etna.core.radiation_damage import radiation_hardness_sweep",
        "radiation_hardness_sweep([8e-4, 10e-4], [5e13, 1e14], [-30.0, -50.0], 1e12)",
    ),
    "dark_current_vs_fluence": (
        "import numpy as np\nfrom etna.core.dark_current import dark_current_vs_fluence",
        "dark_current_vs_fluence(np.geomspace(1e11, 1e13, 6))",
    ),
    "microdosimetric_sweep_devsim_1pt": (
        "from etna.core.optimization import microdosimetric_sweep",
        "microdosimetric_sweep(half_widths_um=[25.0], epi_thicknesses_cm=[10e-4], N_D_bulks=[3.4e13], "
        "V_biases=[20.0], n_lateral_points=5, solver='devsim')",
    ),
    "generation_profiles": (
        "import numpy as np\nfrom etna.core.generation_profiles import alpha_generation_profile, proton_generation_profile, tilted_alpha_generation\nx=np.linspace(0, 30e-4, 2000)",
        "[alpha_generation_profile(x) for _ in range(20)], [proton_generation_profile(x, 62.0) for _ in range(20)]",
    ),
    "lineal_spectrum_1e6": (
        "import numpy as np\nfrom etna.core.microdosimetry import lineal_energy_spectrum\nE=np.random.default_rng(0).lognormal(3,1.2,1_000_000)",
        "lineal_energy_spectrum(E, 20.0)",
    ),
    "compute_kappa_table": (
        "from etna.core.microdosimetry import compute_kappa_table",
        "compute_kappa_table()",
    ),
}


CHILD = r"""
import cProfile, pstats, time, json, sys, io, os, resource
os.environ.setdefault("MPLBACKEND", "Agg")
from etna import set_devsim_quiet
set_devsim_quiet(True)
{setup}
t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
{call}
warm = time.perf_counter() - t0
pr = cProfile.Profile()
t0 = time.perf_counter(); r0 = resource.getrusage(resource.RUSAGE_SELF)
pr.enable()
res = {call}
pr.disable()
wall = time.perf_counter() - t0
r1 = resource.getrusage(resource.RUSAGE_SELF)
cpu = (r1.ru_utime + r1.ru_stime) - (r0.ru_utime + r0.ru_stime)
st = pstats.Stats(pr)
buckets = {{}}
rows = []
for (fn, line, name), (cc, nc, tt, ct, callers) in st.stats.items():
    f = fn.replace("\\", "/")
    if "devsim" in name or "/devsim/" in f:
        b = "devsim"
    elif "/numba/" in f or "Dispatcher" in name or "_compile" in name:
        b = "numba_dispatch"
    elif f.endswith(".py") and ("/etna/core/batched/" in f or "/etna/core/engine2d/" in f) and tt > 0:
        b = "engine_py"
    elif "/numpy/" in f or name.startswith("<method '") and "numpy" in name or "numpy" in name:
        b = "numpy"
    elif "/scipy/" in f:
        b = "scipy"
    elif "/pandas/" in f:
        b = "pandas"
    elif "/matplotlib/" in f:
        b = "matplotlib"
    elif "io.open" in name or "read" in name and fn == "~" or "/csv" in f or "write" in name:
        b = "io"
    elif "/etna/" in f:
        b = "etna_python"
    elif fn == "~":
        b = "builtin_other"
    else:
        b = "other_python"
    buckets[b] = buckets.get(b, 0.0) + tt
    rows.append((tt, ct, nc, f"{{f.split('site-packages/')[-1].split('/etna/')[-1]}}:{{line}}({{name}})"))
rows.sort(reverse=True)
tot = sum(buckets.values()) or 1.0
cum = sorted(((ct, tt, nc, n) for tt, ct, nc, n in rows), reverse=True)
print("__JSON__" + json.dumps(dict(wall=wall, warm=warm, cpu=cpu, profiled_total=tot,
    shares={{k: round(v / tot, 3) for k, v in sorted(buckets.items(), key=lambda kv: -kv[1])}},
    top_tottime=[(round(a, 4), round(b, 4), c, d) for a, b, c, d in rows[:{top}]],
    top_cumtime=[(round(a, 4), round(b, 4), c, d) for a, b, c, d in cum[:{top}]])))
"""


def run(name: str, top: int) -> dict:
    setup, call = PATHS[name]
    code = CHILD.format(setup=setup, call=call, top=top)
    t0 = time.perf_counter()
    p = subprocess.run(
        [sys.executable, "-c", code],
        capture_output=True,
        text=True,
        cwd=HERE.parents[1],
        env={**os.environ, "MPLBACKEND": "Agg"},
    )
    total = time.perf_counter() - t0
    for line in p.stdout.splitlines():
        if line.startswith("__JSON__"):
            d = json.loads(line[8:])
            d["name"] = name
            d["process_total"] = total
            return d
    return {"name": name, "error": p.stderr[-2000:], "process_total": total}


if __name__ == "__main__":
    args = sys.argv[1:]
    top = 15
    if args and args[0] == "--top":
        top = int(args[1])
        args = args[2:]
    names = [n for n in PATHS if not args or any(a in n for a in args)]
    OUT.mkdir(exist_ok=True)
    with open(OUT / "paths.jsonl", "a") as fh:
        for n in names:
            d = run(n, top)
            fh.write(json.dumps(d) + "\n")
            fh.flush()
            if "error" in d:
                print(f"== {n}: ERROR\n{d['error']}")
                continue
            print(
                f"== {n}: wall {d['wall']:.3f}s (warm call {d['warm']:.3f}s) cpu {d['cpu']:.3f}s "
                f"process {d['process_total']:.1f}s shares {d['shares']}"
            )
            for row in d["top_cumtime"][:top]:
                print("   cum", row)
            for row in d["top_tottime"][: max(5, top // 2)]:
                print("   tot", row)

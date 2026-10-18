"""Design optimization utilities for 4H-SiC microdosimeter feasibility study.

Provides:
- Parametric sweep of SV geometry (half-width, epi thickness, doping, bias)
  with a validity gate (full depletion with a 20 % margin + CCE floor) and an
  operating-point ranking (lowest bias first). It does NOT rank by CCE uniformity:
  the modelled device is laterally uniform, so edge/centre CCE is 1 by construction
- Noise floor estimation from dark current (shot noise limit)
- Multi-criteria structure comparison scoring (planar, mesa, 3D, etc.)
- Dark current extraction helper for 2D devices

All units follow project conventions:
    - Lengths: cm (internal), um (user-facing parameters)
    - Currents: A/cm (2D) or A (total)
    - Energies: eV, keV
    - Lineal energy: keV/um

References:
    - Phase 20: charge_collection_2d.py CCE methodology
    - Phase 14: dark_current.py TAT model
    - Phase 23: microdosimetry.py mean_chord_length
    - ICRU Report 36: microdosimetric definitions
"""

import contextlib
import itertools
import logging
import multiprocessing.context as _mpc
import os
import warnings

import devsim
import numpy as np
import pandas as pd

from etna.core.analytical import built_in_potential, full_depletion_voltage_graded
from etna.core.charge_collection_2d import create_2d_dd_device, cce_lateral_scan
from etna.core.dark_current import _resolve_self_consistent, setup_tat_model
from etna.core.devsim_lock import serialized
from etna.core.devsim_reset import reset_devsim_fully
from etna.core.drift_diffusion import extract_steady_state_current
from etna.core.microdosimetry import mean_chord_length
from etna.core.sic_material import SiC4H_Parameters
from etna.core.device import DEFAULT_N_D_NPLUS, DEFAULT_STACK

# AUDIT Mj-3 (v5): minimum CCE for a config to be a trustworthy microdosimeter.
# Below this, edge/center uniformity is meaningless (both-low CCE looks uniform).
CCE_FLOOR = 0.95
# Depletion margin of the validity gate: a row counts as fully depleted only when
# |V_bias| >= V_FD_MARGIN * V_fd. 1.2 covers the +/-20 % epi doping tolerance of
# DESIGN-1 (V_fd scales linearly with N_D at fixed thickness).
V_FD_MARGIN = 1.2

logger = logging.getLogger(__name__)

# Physical constants
Q = 1.602e-19  # C, elementary charge
# E-h pair creation energy: single source of truth in SiC4H_Parameters (7.8 eV,
# measured 4H-SiC W-value). Audit M6: was hardcoded 8.4 eV here, ~7.7% high.
E_PAIR_EV = (
    SiC4H_Parameters().E_pair_eV
)  # eV, electron-hole pair creation energy in 4H-SiC


class _DevsimNonConvergence(RuntimeError):
    """devsim's bias ramp did not converge for this sweep config."""


_OUTSIDE_VALIDITY = (
    "outside model validity: cathode contact on the n- epi (legacy model, no n+ layer), "
    "devsim does not converge here"
)


def _failure_reason(exc):
    """String for the ``failure_reason`` column of a failed sweep row."""
    from etna.core.engine2d import Engine2DOutsideValidity

    if isinstance(exc, Engine2DOutsideValidity):
        return f"[engine2d field threshold] {exc}"
    if isinstance(exc, _DevsimNonConvergence):
        return f"{_OUTSIDE_VALIDITY} [devsim: {exc}]"
    return f"error: {type(exc).__name__}: {exc}"


def analytic_vfd_from_info(params, V_bi):
    """Analytic graded V_fd (V) of a built device's resolved doping (device_info dict).

    Closed form (``full_depletion_voltage_graded``) for the plain exponential profile.
    With the graded n-/n+ transition (``buffer_thickness_cm`` set) the sensitive,
    low-doped epi ends where the transition starts, so V_fd is the moment integral of
    the actual profile (``donor_profile``) from the junction to the TRANSITION START
    (epi end - buffer thickness), via ``full_depletion_voltage_profile``: the same
    convention as the engine's reach-through rule (etna.api.simulation._reach_through_V).
    Depleting the transition itself (N_D up to buffer_end_doping) would take ~130 V on
    the measured-C-V preset and is not "full depletion" of the detector; it is the
    -10..-60 V creep of the measured C-V.
    """
    from etna.core.analytical import full_depletion_voltage_profile
    from etna.core.device import (
        donor_profile,
        donor_table_spec,
        table_sensitive_end_cm,
        tabulated_donor_profile,
    )

    t = params["epi_thickness_cm"]
    if params.get("donor_table") is not None:
        # tabulated profile: moment integral to the table's sensitive-epi end (one decade
        # above its minimum, etna.core.device.table_sensitive_end_cm; a definition, as
        # the transition-start rule below)
        spec = donor_table_spec(t, params["donor_table"])
        y = np.linspace(0.0, table_sensitive_end_cm(spec), 20001)
        return full_depletion_voltage_profile(
            y, tabulated_donor_profile(y, 0.0, spec), V_bi
        )
    if params.get("buffer_thickness_cm") is None:
        return full_depletion_voltage_graded(
            epi_thickness=t,
            N_D_bulk=params["N_D_bulk"],
            N_D_junction=params["N_D_junction"],
            L_transition=params["L_transition"],
            V_bi=V_bi,
        )
    t_low = t - params["buffer_thickness_cm"]
    y = np.linspace(0.0, t_low, 20001)
    nd = donor_profile(
        y,
        0.0,
        params["N_D_junction"],
        params["N_D_bulk"],
        params["L_transition"],
        buffer=(t, params["buffer_thickness_cm"], params["buffer_end_doping"]),
    )
    return full_depletion_voltage_profile(y, nd, V_bi)


@serialized  # lock + quiet devsim log, also inside a spawned worker
def _sweep_one(task):
    """One configuration of ``microdosimetric_sweep`` (module level so it can be pickled)."""
    i, total, hw, epi, nd, vb, n_lateral_points, *rest = task
    solver = rest[0] if rest else "auto"
    # optional extra create_sic_2d_device keywords (n+ layer); absent = builder defaults
    # (the default stack, n+ included)
    dkw = dict(rest[1]) if len(rest) > 1 and rest[1] else {}
    if solver not in ("devsim", "auto"):
        raise ValueError(f"solver must be 'devsim' or 'auto', got {solver!r}")
    logger.info(f"Sweep {i + 1}/{total}: hw={hw}, epi={epi}, nd={nd:.1e}, vb={vb}")
    device_info = None
    used = "devsim"
    failure_reason = ""
    # stack columns of the row (the n+ layer of this config; NaN / False = no n+)
    t_np = dkw.get("nplus_thickness_cm", DEFAULT_STACK["nplus_thickness_cm"])
    stack = {
        "nplus_thickness_um": np.nan if t_np is None else round(t_np * 1e4, 12),
        "N_D_nplus": (
            np.nan if t_np is None else dkw.get("N_D_nplus", DEFAULT_N_D_NPLUS)
        ),
        "nplus_generation": bool(
            t_np is not None and dkw.get("nplus_generation", False)
        ),
    }
    try:
        scan = None
        if solver == "auto":
            # PERF-03 engine: laterally uniform devices inside its validated box only;
            # anything else (or an engine decline) falls back to devsim below.
            from etna.core.engine2d import (
                Engine2DFallback,
                Engine2DOutsideValidity,
                engine2d_unsupported_reason,
                lateral_cce_scan_2d,
            )

            reason = engine2d_unsupported_reason(
                hw, vb, epi, nd, n_lateral_points, **dkw
            )
            if reason is None:
                try:
                    used = "engine2d"
                    scan = lateral_cce_scan_2d(
                        hw,
                        vb,
                        epi_thickness_cm=epi,
                        N_D_bulk=nd,
                        n_lateral_points=n_lateral_points,
                        **dkw,
                    )
                    params = scan["device_info"]
                except Engine2DOutsideValidity:
                    # The same failed (NaN) row devsim produces there (its ramp does not
                    # converge): outside model validity (cathode on the n- epi), not physics.
                    raise
                except Engine2DFallback as e:
                    used = "devsim"
                    logger.info(
                        f"Sweep {i + 1}/{total}: engine2d declined ({e}); using devsim"
                    )
                except (
                    Exception
                ) as e:  # any other engine error must not become a quiet NaN row
                    used = "devsim"
                    scan = None
                    logger.warning(
                        f"Sweep {i + 1}/{total}: engine2d error ({e!r}); using devsim"
                    )
            else:
                logger.info(
                    f"Sweep {i + 1}/{total}: engine2d unsupported ({reason}); using devsim"
                )
        if scan is None:
            try:
                device_info = create_2d_dd_device(
                    half_width_um=hw,
                    V_bias=vb,
                    epi_thickness_cm=epi,
                    N_D_bulk=nd,
                    **dkw,
                )
            except RuntimeError as e:
                if "ramp_bias: failed to converge" in str(e):
                    if t_np is not None:
                        raise  # n+ device: no contact-on-epi wall, a plain failure
                    raise _DevsimNonConvergence(str(e)) from e
                raise
            scan = cce_lateral_scan(device_info, n_points=n_lateral_points)
            params = device_info
        cce_vals = np.array(scan["cce_values"])

        # AUDIT Mj-3: full-depletion gate. A partially-depleted device with
        # both-low CCE shows a deceptively uniform edge/center ratio and
        # would rank artificially high. Compute V_fd from the device's
        # actual graded profile (NOT the uniform N_D_bulk estimate, which is
        # anti-conservative) and require both full depletion AND an absolute
        # CCE floor before the uniformity ratio is trustworthy.
        V_bi = built_in_potential(
            params["N_A_ionized"],
            params["N_D_junction"],
            params["n_i"],
        )
        V_fd = analytic_vfd_from_info(params, V_bi)
        raw_center_cce = cce_vals[0]
        raw_edge_cce = cce_vals[-1]
        # Depletion margin (V_FD_MARGIN = 1.2): covers the +/-20 % doping tolerance.
        is_fully_depleted = abs(vb) >= V_FD_MARGIN * V_fd
        if raw_center_cce is None or raw_edge_cce is None:
            # Defensive: a degenerate/undepleted config can make
            # cce_lateral_scan() yield a non-numeric (None) CCE entry
            # instead of a float or NaN. Treat that explicitly as a
            # floor failure rather than letting float() below raise a
            # TypeError that gets swallowed by the generic except block.
            center_cce = np.nan
            edge_cce = np.nan
            passes_cce_floor = False
        else:
            center_cce = float(raw_center_cce)
            edge_cce = float(raw_edge_cce)
            passes_cce_floor = center_cce >= CCE_FLOOR and edge_cce >= CCE_FLOOR
        record = {
            "half_width_um": hw,
            "epi_thickness_cm": epi,
            "N_D_bulk": nd,
            "V_bias": vb,
            "center_cce": center_cce,
            "edge_cce": edge_cce,
            "edge_center_ratio": float(scan["edge_to_center_ratio"]),
            "cce_std": float(np.std(cce_vals)),
            "V_fd": float(V_fd),
            "is_fully_depleted": bool(is_fully_depleted),
            "passes_cce_floor": bool(passes_cce_floor),
            "is_valid": bool(is_fully_depleted and passes_cce_floor),
            "solver": used,
            "failure_reason": "",
            **stack,
        }
    except Exception as e:
        failure_reason = _failure_reason(e)
        logger.warning(
            f"Sweep {i + 1}/{total} failed (hw={hw}, epi={epi}, "
            f"nd={nd:.1e}, vb={vb}): {failure_reason}"
        )
        record = {
            "half_width_um": hw,
            "epi_thickness_cm": epi,
            "N_D_bulk": nd,
            "V_bias": vb,
            "center_cce": np.nan,
            "edge_cce": np.nan,
            "edge_center_ratio": np.nan,
            "cce_std": np.nan,
            "V_fd": np.nan,
            "is_fully_depleted": False,
            "passes_cce_floor": False,
            "is_valid": False,
            "solver": used,
            "failure_reason": failure_reason,
            **stack,
        }
        # Reset devsim global state after failure to prevent
        # contamination of subsequent configurations (a failed
        # ramp_bias can leave solver state that causes all later
        # equilibrium solves to diverge). reset_devsim_fully also clears
        # any cylindrical-axis globals (P03/P30) and preserves the solver.
        reset_devsim_fully(preserve_solver=True)
    finally:
        if device_info is not None:
            try:
                devsim.delete_device(device=device_info["device_name"])
            except Exception:
                pass
    return record


# Thread-pool sizes of the native libraries in a sweep worker: one thread each, so
# n_jobs workers use n_jobs cores (numba's default is one thread per core per worker).
_WORKER_THREAD_ENV = (
    "NUMBA_NUM_THREADS",
    "OMP_NUM_THREADS",
    "OPENBLAS_NUM_THREADS",
    "MKL_NUM_THREADS",
    "VECLIB_MAXIMUM_THREADS",
)


def _worker_init():
    """Spawn-pool initializer: one thread per native pool in the worker.

    The environment caps are already inherited from the parent (``_spawn_pool``), which
    is what counts for OpenBLAS / OpenMP: they size their pools when the worker first
    imports numpy / numba, i.e. while unpickling this initializer, before it runs. They
    are repeated here for any process the worker starts, and numba's runtime count is
    set explicitly.
    """
    for key in _WORKER_THREAD_ENV:
        os.environ[key] = "1"
    import numba

    numba.set_num_threads(1)


@contextlib.contextmanager
def _main_module_hidden():
    """Show the spawn machinery a ``__main__`` without ``__file__`` / ``__spec__``.

    A spawned child re-executes the parent's ``__main__`` file (as ``__mp_main__``)
    unless the main module has neither. Under streamlit (the app, and AppTest in the
    tests) ``__main__`` is the page script, installed by the script runner and never
    restored: the child re-ran it and died (BrokenProcessPool; the AppTest temp
    script raised NameError). The sweep workers only run importable etna functions.
    The swap lasts one process launch (a few ms) and the parent's module is put back.
    """
    import sys
    import types

    saved = sys.modules.get("__main__")
    sys.modules["__main__"] = types.ModuleType("__main__")
    try:
        yield
    finally:
        sys.modules["__main__"] = saved


class _NoMainSpawnProcess(_mpc.SpawnProcess):
    """SpawnProcess launched without the parent's ``__main__`` (module level: the
    process object is pickled to the child)."""

    @staticmethod
    def _Popen(process_obj):
        with _main_module_hidden():
            return _mpc.SpawnProcess._Popen(process_obj)


class _NoMainSpawnContext(_mpc.SpawnContext):
    """The "spawn" context, its workers launched without the parent's ``__main__``."""

    Process = _NoMainSpawnProcess


def _no_main_spawn_context():
    return _NoMainSpawnContext()


@contextlib.contextmanager
def _spawn_pool(max_workers):
    """ProcessPoolExecutor (spawn) whose workers run single-threaded native pools.

    The caps are put in ``os.environ`` while the workers are started (a spawned child
    inherits the parent's environment) and restored afterwards; the parent's own pools
    are already sized, so it is not affected. The caller holds the devsim lock
    (``serialized``), so no other etna call starts processes meanwhile. The workers do
    not import the parent's ``__main__`` (``_main_module_hidden``): submitted functions
    must live in importable modules.
    """
    from concurrent.futures import ProcessPoolExecutor

    saved = {key: os.environ.get(key) for key in _WORKER_THREAD_ENV}
    os.environ.update({key: "1" for key in _WORKER_THREAD_ENV})
    try:
        with ProcessPoolExecutor(
            max_workers=max_workers,
            mp_context=_no_main_spawn_context(),
            initializer=_worker_init,
        ) as pool:
            yield pool
    finally:
        for key, value in saved.items():
            if value is None:
                os.environ.pop(key, None)
            else:
                os.environ[key] = value


@serialized
def microdosimetric_sweep(
    half_widths_um=(25, 50, 100, 150),
    epi_thicknesses_cm=(5e-4, 10e-4, 20e-4),
    N_D_bulks=(5e13, 8.5e13, 5e14),
    V_biases=(20, 40, 60),
    n_lateral_points=10,
    max_configs=None,
    n_jobs=1,
    solver="auto",
    nplus_thickness_cm=DEFAULT_STACK["nplus_thickness_cm"],
    N_D_nplus=None,
    nplus_generation=False,
    substrate_thickness_cm=DEFAULT_STACK["substrate_thickness_cm"],
    N_D_junction=None,
    L_transition=None,
    buffer_thickness_cm=DEFAULT_STACK["buffer_thickness_cm"],
    buffer_end_doping=DEFAULT_STACK["buffer_end_doping"],
    donor_table=None,
):
    """Sweep SV geometry and doping parameters with a validity gate and operating-point ranking.

    Creates a full combinatorial grid of (half_width, epi_thickness, N_D_bulk,
    V_bias) and for each configuration runs a lateral CCE scan. The modelled
    device is laterally uniform, so edge_center_ratio is 1 up to solver noise and
    is NOT a uniformity result or a ranking key (see _rank_sweep_results). CCE is
    relative to the generation inside the simulated epi (in-epi generation), not
    to the full alpha deposit.

    Parameters
    ----------
    half_widths_um : tuple of float
        SV half-widths to sweep (micrometers).
    epi_thicknesses_cm : tuple of float
        Epitaxial layer thicknesses to sweep (cm).
    N_D_bulks : tuple of float
        Bulk donor concentrations to sweep (cm^-3).
    V_biases : tuple of float
        Reverse bias voltages to sweep (V).
    n_lateral_points : int
        Number of lateral positions per CCE scan.
    max_configs : int or None
        If set, truncate the parameter grid to this many configurations
        (with a warning).
    n_jobs : int
        Worker processes (default 1 = serial in this process). Configurations are
        independent fresh builds, so results are identical to the serial run.
    solver : {"auto", "devsim"}
        "auto" (default) uses the numba
        2D engine (etna.core.engine2d, PERF-03) for configs it supports -- laterally
        uniform default devices inside its validated box -- and devsim for the rest
        or whenever the engine declines. "devsim" builds every config with devsim.
        Configs outside model validity (the cathode is an ohmic contact on the n- epi
        with no n+ substrate; past reach-through devsim's ramp stops converging at a
        fixed contact-cell field, which is not physical punch-through) are reported as
        the same failed (NaN) row, with ``failure_reason`` set.
    nplus_thickness_cm, N_D_nplus, nplus_generation :
        n+ layer under the epi for every config (create_sic_2d_device keywords).
        Default: the default stack (1 um, N_D_nplus None = 1e19, epi-only generation:
        matches the measurement geometry, alphas tilted 65 deg with projected depth
        ~7.6 um never reach the n+). nplus_thickness_cm=None = the LEGACY device
        (cathode on the n- epi, non-convergence wall; pass substrate_thickness_cm=1e-4
        for the legacy p+). With n+ the cathode sits on
        the n+ and there is no contact-on-epi wall. "auto" serves n+ configs inside
        etna.core.engine2d.ENGINE2D_NPLUS_BOX (parity vs devsim:
        lab/solver2d/notes/LOG.md "2D engine with n+"); outside it, devsim.
    substrate_thickness_cm, N_D_junction, L_transition :
        p+ anode thickness (cm; default the stack's 0.3 um) and graded epi doping of
        every config (None = create_sic_2d_device default = DEFAULT_STACK doping),
        forwarded with the n+ keywords.
    buffer_thickness_cm, buffer_end_doping :
        Graded n-/n+ transition of every config (create_sic_2d_device keywords; default
        DEFAULT_STACK, both None = abrupt epi end). The V_fd gate then integrates the
        actual profile to the transition start (``analytic_vfd_from_info``). "auto"
        serves it inside ENGINE2D_BUFFER_BOX. The default transition (Lb 1.14 um,
        N_end 4.5e16) was FITTED to the C-V of the 10 um epi; it is applied unchanged
        (same Lb before the epi end, same end doping) to every epi thickness of the
        grid, so the rows at other thicknesses (e.g. the default 5 and 20 um) are a
        model extrapolation, not a calibrated device.
    donor_table : tuple or None
        Tabulated epi profile ((depth_um, N_D), ...; create_sic_2d_device keyword) for
        every config: it REPLACES the graded doping, so the ``N_D_bulks`` grid no longer
        changes the device (the column is kept as the label of the row) and the V_fd gate
        integrates the table to its sensitive-epi end (``analytic_vfd_from_info``). The
        table must cover every epi thickness of the grid. "auto" serves it inside
        ENGINE2D_TABLE_BOX.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: half_width_um, epi_thickness_cm, N_D_bulk,
        V_bias, center_cce, edge_cce, edge_center_ratio, cce_std, solver
        ("devsim" or "engine2d": which solver produced the row), plus the
        Mj-3 validity gate: V_fd (full-depletion voltage from the graded
        profile), is_fully_depleted (|V_bias| >= 1.2 * V_fd: 20 % margin for the
        +/-20 % doping tolerance of DESIGN-1), passes_cce_floor (both center &
        edge CCE >= CCE_FLOOR=0.95), and is_valid (both); and failure_reason
        (str: empty for a solved row; otherwise outside model validity -- engine
        field threshold or devsim non-convergence -- or the exception type and
        message). Ranked by _rank_sweep_results: valid rows first, then lowest
        V_bias, then largest V_bias - V_fd; NOT by edge_center_ratio. Invalid
        (partially depleted, low-CCE or failed) configs are RETAINED and
        flagged, not dropped.
    """
    grid = list(
        itertools.product(half_widths_um, epi_thicknesses_cm, N_D_bulks, V_biases)
    )
    total = len(grid)

    if max_configs is not None and total > max_configs:
        warnings.warn(
            f"Parameter grid has {total} configurations, truncating to {max_configs}. "
            f"Reduce sweep ranges for exhaustive coverage.",
            UserWarning,
            stacklevel=2,
        )
        grid = grid[:max_configs]
        total = len(grid)

    # nplus_thickness_cm is always passed: None must reach the builder as an explicit
    # "no n+" (its default is the stack's n+ layer).
    dkw = {"nplus_thickness_cm": nplus_thickness_cm}
    if nplus_thickness_cm is not None:
        dkw["nplus_generation"] = bool(nplus_generation)
        if N_D_nplus is not None:
            dkw["N_D_nplus"] = N_D_nplus
    # p+ anode / graded doping of every config (None = create_sic_2d_device default)
    for key, value in (
        ("substrate_thickness_cm", substrate_thickness_cm),
        ("N_D_junction", N_D_junction),
        ("L_transition", L_transition),
    ):
        if value is not None:
            dkw[key] = value
    # the transition keys are always passed: None must reach the builder as an explicit
    # "abrupt epi end" (its default is the stack's graded transition)
    dkw["buffer_thickness_cm"] = buffer_thickness_cm
    dkw["buffer_end_doping"] = buffer_end_doping
    if donor_table is not None:
        dkw["donor_table"] = tuple((float(d), float(n)) for d, n in donor_table)
    tasks = [
        (i, total, hw, epi, nd, vb, n_lateral_points, solver) + ((dkw,) if dkw else ())
        for i, (hw, epi, nd, vb) in enumerate(grid)
    ]
    if solver not in ("devsim", "auto"):
        raise ValueError(f"solver must be 'devsim' or 'auto', got {solver!r}")
    if n_jobs is not None and n_jobs > 1 and total > 1:
        # devsim state is per process and every config is a fresh build, so the
        # results are identical to the serial loop (executor.map keeps the order).
        with _spawn_pool(min(n_jobs, total)) as pool:
            records = list(pool.map(_sweep_one, tasks))
    else:
        records = [_sweep_one(t) for t in tasks]

    return _rank_sweep_results(pd.DataFrame(records))


def _rank_sweep_results(df):
    """Rank microdosimetric-sweep results with the Mj-3 validity gate.

    Valid configs (fully-depleted AND above the CCE floor) rank above all invalid
    ones, which are retained (flagged is_valid=False), not dropped. Within the valid
    group the operating point is ranked: lowest V_bias first (most margin from the
    model's validity wall -- devsim non-convergence with the cathode on the n- epi, not a
    physical limit -- and there is no breakdown model, so a lower field is the safe side),
    then the largest full-depletion margin V_bias - V_fd, then half_width_um ascending
    as a deterministic tie-break only (it is NOT physics).

    edge_center_ratio is deliberately NOT a ranking key any more: the modelled device
    is laterally uniform (depth-only doping, full-width contacts, mirror sides), so the
    ratio is 1 up to solver noise (~1e-6) and used to order that noise (PERF-03). The
    column is kept as a sanity check: it should stay within ~1e-5 of 1. Extra columns
    (e.g. failure_reason, a string) are carried through untouched. Pure function.
    """
    df = df.copy()
    keys, ascending = ["is_valid"], [False]
    if "V_bias" in df.columns:
        keys.append("V_bias")
        ascending.append(True)
        if "V_fd" in df.columns:
            df["_margin"] = df["V_bias"] - df["V_fd"]
            keys.append("_margin")
            ascending.append(False)
    if "half_width_um" in df.columns:
        keys.append("half_width_um")
        ascending.append(True)
    df = df.sort_values(keys, ascending=ascending, na_position="last", kind="stable")
    return df.drop(columns=["_margin"], errors="ignore").reset_index(drop=True)


def estimate_noise_floor(
    I_dark_A,
    t_shaping_s=1e-6,
    sv_thickness_um=10.0,
    sv_width_um=100.0,
    sv_depth_um=None,
    k_sigma=3.0,
):
    """Estimate detector-intrinsic noise floor from dark current.

    Computes the minimum detectable charge, energy, and lineal energy
    from shot noise on the dark current. This represents the fundamental
    detector-intrinsic limit; real systems will have higher noise from
    readout electronics, electromagnetic interference, and other sources.

    Parameters
    ----------
    I_dark_A : float
        Dark current in Amperes (total, not density).
    t_shaping_s : float
        Shaping time of the readout electronics (seconds).
    sv_thickness_um : float
        Sensitive volume thickness (micrometers).
    sv_width_um : float
        Sensitive volume width (micrometers).
    sv_depth_um : float or None
        Sensitive volume depth (micrometers). If None, uses slab
        approximation for mean chord length.
    k_sigma : float
        Number of standard deviations for detection threshold
        (default 3 for 99.7% confidence).

    Returns
    -------
    dict
        Dictionary with keys:
        - I_dark_A: input dark current
        - sigma_shot_C: shot noise RMS charge (Coulombs)
        - Q_min_fC: minimum detectable charge (femtocoulombs)
        - E_min_keV: minimum detectable energy (keV)
        - y_min_keV_um: minimum detectable lineal energy (keV/um)
        - l_bar_um: mean chord length (um)
        - t_shaping_s: shaping time used
        - k_sigma: detection threshold used
    """
    sigma_shot = np.sqrt(2.0 * Q * abs(I_dark_A) * t_shaping_s)
    Q_min = k_sigma * sigma_shot

    # Energy threshold from pair creation energy
    E_min_eV = Q_min * E_PAIR_EV / Q
    E_min_keV = E_min_eV / 1e3

    # Mean chord length for lineal energy conversion
    l_bar = mean_chord_length(sv_thickness_um, sv_width_um, sv_depth_um)

    # Minimum detectable lineal energy
    y_min_keV_um = E_min_keV / l_bar

    return {
        "I_dark_A": I_dark_A,
        "sigma_shot_C": float(sigma_shot),
        "Q_min_fC": float(Q_min * 1e15),
        "E_min_keV": float(E_min_keV),
        "y_min_keV_um": float(y_min_keV_um),
        "l_bar_um": float(l_bar),
        "t_shaping_s": t_shaping_s,
        "k_sigma": k_sigma,
    }


def score_structures(
    metrics_dict,
    weights=None,
):
    """Score and rank alternative microdosimeter structures.

    Applies multi-criteria weighted scoring with min-max normalization
    across structures. Metrics are classified as "higher is better"
    (cce_uniformity) or "lower is better" (noise_floor, spectral_resolution,
    fabrication_complexity).

    Parameters
    ----------
    metrics_dict : dict
        {structure_name: {metric_name: value}} where metrics are:
        - cce_uniformity: edge/center CCE ratio (higher = better)
        - noise_floor: y_min in keV/um (lower = better)
        - spectral_resolution: y_D/y_F ratio or d(y) peak std
          (lower = better)
        - fabrication_complexity: 1-4 scale (lower = better)
    weights : dict or None
        {metric_name: weight}. Default: cce_uniformity=0.30,
        noise_floor=0.20, spectral_resolution=0.20,
        fabrication_complexity=0.30. Weights must sum to 1.0.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns: structure, each raw metric, each
        normalized metric (suffixed _norm), weighted_score.
        Sorted by weighted_score descending (best first).
    """
    if weights is None:
        weights = {
            "cce_uniformity": 0.30,
            "noise_floor": 0.20,
            "spectral_resolution": 0.20,
            "fabrication_complexity": 0.30,
        }

    # Metrics where higher is better
    higher_is_better = {"cce_uniformity"}

    # Build raw DataFrame
    structures = list(metrics_dict.keys())
    metric_names = list(weights.keys())

    rows = []
    for name in structures:
        row = {"structure": name}
        for m in metric_names:
            row[m] = metrics_dict[name].get(m, 0.0)
        rows.append(row)

    df = pd.DataFrame(rows)

    # Normalize each metric to [0, 1]
    for m in metric_names:
        vals = df[m].values.astype(float)
        v_min = vals.min()
        v_max = vals.max()

        if v_max == v_min:
            # All equal -- assign 0.5
            df[f"{m}_norm"] = 0.5
        else:
            normalized = (vals - v_min) / (v_max - v_min)
            if m not in higher_is_better:
                # Lower is better: invert
                normalized = 1.0 - normalized
            df[f"{m}_norm"] = normalized

    # Compute weighted score
    df["weighted_score"] = sum(weights[m] * df[f"{m}_norm"] for m in metric_names)

    df = df.sort_values("weighted_score", ascending=False).reset_index(drop=True)
    return df


@serialized
def get_dark_current_2d(half_width_um=50.0, V_bias=50.0, **device_kwargs):
    """Extract total dark current from a 2D DD device with TAT model.

    Creates a 2D drift-diffusion device, applies the trap-assisted tunneling
    model, re-solves, and extracts the dark current from the cathode contact.
    The 2D current (A/cm) is multiplied by the SV depth to get total current
    in Amperes.

    Parameters
    ----------
    half_width_um : float
        SV half-width (micrometers).
    V_bias : float
        Reverse bias voltage (V, positive).
    **device_kwargs
        Additional arguments for create_2d_dd_device (e.g.,
        epi_thickness_cm, N_D_bulk).

    Returns
    -------
    float
        Dark current in Amperes.
    """
    device_info = None
    try:
        device_info = create_2d_dd_device(
            half_width_um=half_width_um,
            V_bias=V_bias,
            **device_kwargs,
        )

        # Apply the TAT model and solve to the self-consistent frozen-field state (as the
        # 1D dark_current_sweep; vs the former single solve: 7e-15 at hw 25, 50 V, +0.3 s)
        setup_tat_model(device_info)
        _resolve_self_consistent(device_info["device_name"], device_info["region_name"])

        # Extract 2D current (A/cm)
        I_2d = extract_steady_state_current(device_info, contact="cathode")

        # Convert to total current. The 2D domain is the HALF device (x from 0 to
        # half_width, mirror boundaries), so the full SV has twice that current per cm
        # of depth; then multiply by the SV depth (square cross-section, 2 * half_width).
        # Missing the factor 2 made this 0.501x the 1D-twin J times the SV area (PERF-03).
        sv_depth_cm = 2.0 * half_width_um * 1e-4
        I_total = 2.0 * I_2d * sv_depth_cm

        logger.info(
            f"Dark current 2D: I_2d={I_2d:.3e} A/cm, "
            f"sv_depth={sv_depth_cm * 1e4:.0f} um, I_total={I_total:.3e} A"
        )

        return float(I_total)

    finally:
        if device_info is not None:
            try:
                devsim.delete_device(device=device_info["device_name"])
            except Exception:
                pass

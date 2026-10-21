"""Facades of the single-ion track engine (``etna.core.track``, issue #9).

* ``run_track_event``: one ion track at normal incidence (axisymmetric transient):
  cathode current vs time, CCE, recombined fraction, ledgers;
* ``run_track_cce_let``: CCE / pulse-height defect vs LET on one device and bias;
* ``run_tilted_alpha_track``: CCE vs reverse bias of a tilted alpha with the plasma
  (columnar) recombination of the track (2.5D: 1D column linear response + transverse
  plasma slices along the path, ``etna.core.track.tilted``).

All three are opt-in physics: no other facade calls them, and no default changes. The
track engine is numba (no devsim). Results are ``SimResult`` envelopes.
"""

from __future__ import annotations

import numpy as np

from etna.api.device import DeviceConfig
from etna.api.results import SimResult

MESHES = ("tiny", "fast", "default")


def _radial(radial, s_nm, beta=None):
    from etna.core.track.radial import RadialModel, beta_of, gaussian, thermalized, whk

    if isinstance(radial, RadialModel) or callable(radial):
        return radial
    if radial == "gaussian":
        return gaussian(s_nm * 1e-7)
    if radial == "whk":
        b = beta_of(1.0) if beta is None else beta
        return thermalized(whk(b), s_nm * 1e-7) if s_nm > 0 else whk(b)
    raise ValueError("radial must be 'gaussian', 'whk' or a RadialModel")


def _device(config, mesh, auger_model, carrier_scattering, auger=True):
    from etna.core.track import phd
    from etna.core.track.grid import device_from_config, radial_nodes

    if mesh not in MESHES:
        raise ValueError(f"mesh must be one of {MESHES}")
    m = dict({"tiny": phd.MESH_TINY, "fast": phd.MESH_FAST,
              "default": phd.MESH_DEFAULT}[mesh])  # fmt: skip
    r = radial_nodes(m.pop("r_min_cm"), m.pop("R_cm"), m.pop("r_growth"))
    return device_from_config(
        config or DeviceConfig(), r=r, auger=auger, auger_model=auger_model,
        carrier_scattering=carrier_scattering, **m,
    )  # fmt: skip


def run_track_event(
    config: DeviceConfig = None,
    LET_keV_um: float = 100.0,
    V_bias: float = 50.0,
    *,
    radial="whk",
    s_nm: float = 10.0,
    mesh: str = "fast",
    region: str = "epi",
    auger_model: str = "constant",
    carrier_scattering: bool = False,
    rtol: float = 1e-4,
    excess=None,
) -> SimResult:
    """One ion at normal incidence: x = time (s), y = cathode signal current (e/s).

    ``radial``: "whk" (WHK at 1 MeV/u thermalized by ``s_nm``; s_nm = 0: bare),
    "gaussian" (width ``s_nm``) or a ``RadialModel``. ``V_bias``: reverse bias
    magnitude (V). metadata: cce (pair balance), recombined_fraction, wrong_contact,
    ledger_residual, ledger_floor, converged, n_steps, wall_s, excess (unknowns form),
    peak_excess_cm3. The thermalization radius of 4H-SiC is NOT known: the plasma loss
    is conditional on ``radial`` (docs/USAGE.md section 19)."""
    from etna.core.track.phd import simulate_track

    dev = _device(config, mesh, auger_model, carrier_scattering)
    r = simulate_track(device=dev, LET_keV_um=LET_keV_um, V_bias=abs(V_bias),
                       radial=_radial(radial, s_nm), region=region, rtol=rtol,
                       excess=excess)  # fmt: skip
    md = dict(
        cce=r.cce,
        recombined_fraction=r.recombined_fraction,
        wrong_contact=r.wrong_contact / r.N_gen,
        N_gen=r.N_gen,
        ledger_residual=r.ledger_residual,
        ledger_floor=r.extra["ledger_floor"],
        converged=r.converged,
        n_steps=r.n_steps,
        wall_s=r.wall_s,
        excess=r.extra.get("excess"),
        peak_excess_cm3=r.peak_excess,
        LET_keV_um=LET_keV_um,
        V_bias=abs(V_bias),
        mesh=mesh,
        solver="track (axisymmetric numba DD)",
        label="predicted (conditional on the radial track model)",
    )
    return SimResult(
        config or DeviceConfig(), "track_event", r.times, r.current_cathode, md
    )


def run_track_cce_let(
    config: DeviceConfig = None,
    LET_values=(1.0, 10.0, 100.0, 1000.0),
    V_bias: float = 50.0,
    *,
    radial="gaussian",
    s_nm: float = 50.0,
    mesh: str = "fast",
    region: str = "epi",
    auger_model: str = "constant",
    carrier_scattering: bool = False,
    rtol: float = 1e-4,
) -> SimResult:
    """CCE vs LET (keV/um) of normal-incidence tracks on one device / dark state:
    x = LET, y = CCE; metadata ``table`` (the ``cce_vs_let`` DataFrame as a dict of
    lists), PHD, recombined fraction. Default radial: Gaussian 50 nm (commercial style).
    """
    from etna.core.track.phd import simulate_track
    from etna.core.track.solver import dark_state

    dev = _device(config, mesh, auger_model, carrier_scattering)
    dark = dark_state(dev, abs(V_bias))
    rows = []
    for let in LET_values:
        r = simulate_track(device=dev, dark=dark, LET_keV_um=float(let),
                           V_bias=abs(V_bias), radial=_radial(radial, s_nm),
                           region=region, rtol=rtol, keep_series=False)  # fmt: skip
        rows.append((float(let), r.cce, r.recombined_fraction, r.ledger_residual,
                     r.converged, r.wall_s))  # fmt: skip
    a = np.array(rows, dtype=float)
    md = dict(
        PHD=1.0 - a[:, 1],
        recombined_fraction=a[:, 2],
        ledger_residual=a[:, 3],
        converged=a[:, 4].astype(bool),
        wall_s=a[:, 5],
        V_bias=abs(V_bias),
        mesh=mesh,
        solver="track (axisymmetric numba DD)",
        label="predicted (conditional on the radial track model)",
    )
    return SimResult(config or DeviceConfig(), "track_cce_let", a[:, 0], a[:, 1], md)


def run_tilted_alpha_track(
    config: DeviceConfig = None,
    V_values=(0.0, 5.0, 10.0, 20.0, 30.0, 40.0, 50.0),
    *,
    theta_deg: float = 65.0,
    alpha_range_um: float = 18.0,
    E0_MeV: float = 5.5,
    dead_layer_um: float = 0.0,
    tau_scale: float = 1.0,
    s_nm: float = 10.0,
    R_um: float = 10.0,
    auger_model: str = "constant",
    carrier_scattering: bool = False,
    mu_field: bool = False,
    workers: int = 1,
    table=None,
) -> SimResult:
    """CCE vs reverse bias of a tilted alpha WITH track-plasma recombination (2.5D).

    x = -V (the run_cce_tilted_alpha sign convention), y = CCE (relative to the epi
    pairs). metadata: cce_column (1D linear response), plasma_loss, shape (normalised
    at the largest |V|), table (the PlasmaTable as a dict), notes (approximations).

    ``tau_scale`` multiplies both SRH lifetimes (column and plasma slices). With the
    library lifetimes (tau_scale 1) Fig. 4b is NOT reproduced; the -20..-50 V tail
    needs tau_scale ~1/20-1/30 (a CALIBRATION, lab/track_plasma/notes/LOG.md, not a
    preset: it misses 0..-10 V and is ~100x below good-epi ambipolar lifetimes).

    ``table``: a precomputed ``PlasmaTable`` (its physics must match the arguments).
    Cost: one plasma table (~80 transverse slices, ~1-3 s each; ``workers`` > 1 runs
    them in processes) plus a column transient per bias."""
    from etna.api.device import microdosimeter_config
    from etna.core.track.tilted import tilted_alpha_cce

    cfg = config if config is not None else microdosimeter_config()
    r = tilted_alpha_cce(
        cfg, V_values, theta_deg=theta_deg, alpha_range_um=alpha_range_um,
        E0_MeV=E0_MeV, dead_layer_um=dead_layer_um, tau_scale=tau_scale,
        s_nm=s_nm, workers=workers, table=table,
        table_kw=dict(R_um=R_um, auger_model=auger_model,
                      carrier_scattering=carrier_scattering, mu_field=mu_field),
    )  # fmt: skip
    md = dict(
        cce_column=r.cce_column,
        plasma_loss=r.plasma_loss,
        shape=r.shape(float(np.max(r.V))),
        table=r.table.to_dict(),
        notes=r.notes,
        column_ledger=r.column_ledger,
        tau_scale=tau_scale,
        s_nm=s_nm,
        solver="track 2.5D (column + transverse plasma slices, numba DD)",
        label="predicted (conditional on tau_scale and the radial model)",
    )
    return SimResult(cfg, "cce_tilted_track", -r.V, r.cce, md)

__all__ = [
    "run_tilted_alpha_track",
    "run_track_cce_let",
    "run_track_event",
]

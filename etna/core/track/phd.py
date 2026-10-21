"""Single-ion track: generation along the path, CCE / pulse-height defect vs LET.

High-level entry points of the axisymmetric track engine:

* ``track_pairs``: pairs per finite volume for a straight track on the axis, line
  density LET / E_pair along z, radial shape from ``etna.core.track.radial`` (lumped
  per annulus exactly);
* ``simulate_track``: one event (device + dark state + transient), a TrackResult;
* ``cce_vs_let``: CCE, recombined fraction and pulse-height defect over a LET grid on
  one device / dark state (the plasma-regime counterpart of
  ``etna.core.single_particle.build_cce_let_table``, whose 2D Cartesian sheet track
  stays linear by construction);
* ``ion_let_profile``: LET(depth) and E/u(depth) of a proton or alpha from the NIST
  PSTAR / ASTAR tables (``etna.core.stopping_power``) for a Bragg-curve track.

Normal incidence only (the track is the symmetry axis). The pulse-height defect here
is the RECOMBINATION part only (plasma / columnar recombination by SRH + Auger in the
dense track) plus any charge not collected; window / dead-layer and nuclear-stopping
defects are not modelled.
"""

from __future__ import annotations

import time

import numpy as np
import pandas as pd

from etna.core.sic_material import SiC4H_Parameters

from .grid import device_from_config, radial_nodes
from .radial import beta_of, gaussian, whk
from .solver import TrackSystem, dark_state, run_transient

E_PAIR_EV = SiC4H_Parameters().E_pair_eV

# Mesh presets (kwargs of device_from_config / radial_nodes). FAST: tests and quick
# scans (21 x 110 nodes); DEFAULT: the bench numbers (~31 x 150).
MESH_FAST = dict(r_min_cm=1e-7, R_cm=20e-4, r_growth=1.6, h_min_cm=5e-7,
                 h_max_cm=0.5e-4, z_growth=1.4)  # fmt: skip
MESH_TINY = dict(
    r_min_cm=2e-7, R_cm=10e-4, r_growth=2.0, h_min_cm=1e-6, h_max_cm=1e-4, z_growth=1.6
)  # fmt: skip  (unit tests only)
MESH_DEFAULT = dict(r_min_cm=0.5e-7, R_cm=30e-4, r_growth=1.4, h_min_cm=3e-7,
                    h_max_cm=0.3e-4, z_growth=1.3)  # fmt: skip


def build_track_device(config, mesh=None, auger=True):
    """TrackDevice of ``config`` on the mesh preset / dict ``mesh`` (None = DEFAULT)."""
    m = dict(MESH_DEFAULT if mesh is None else mesh)
    r = radial_nodes(m.pop("r_min_cm"), m.pop("R_cm"), m.pop("r_growth"))
    return device_from_config(config, r=r, auger=auger, **m)


def track_pairs(dev, LET_keV_um, radial=None, region="epi", E_pair_eV=E_PAIR_EV):
    """(N,) pairs per node of a track on the axis entering at the anode (z = 0).

    ``LET_keV_um``: a number, or a callable of the depth z (cm from the anode surface)
    returning keV/um (0 beyond the range). ``radial``: a RadialModel, a callable of z
    returning one (speed-dependent track structure), or None = gaussian(50 nm).
    ``region``: "epi" (pairs only in the n- epi, the build_cce_let_table convention:
    the p+ is part of the dead layer) or "all" (the whole path, p+ and n+ included).
    Contact nodes get none. Each row j receives LET(z_j)/E_pair * dz_j pairs (dz_j the
    axial dual length), split over the annuli by the radial CDF.
    """
    if region not in ("epi", "all"):
        raise ValueError("region must be 'epi' or 'all'")
    radial = gaussian() if radial is None else radial
    z = dev.z
    let = LET_keV_um(z) if callable(LET_keV_um) else np.full(z.size, float(LET_keV_um))
    let = np.asarray(let, dtype=float)
    line = let * 1e7 / E_pair_eV  # pairs/cm
    if region == "epi":
        line = np.where((z > dev.junction_cm) & (z < dev.epi_end_cm), line, 0.0)
    line[0] = line[-1] = 0.0  # contact rows
    per_row = line * dev.axial_dual()
    faces = dev.radial_faces()
    pairs = np.zeros((dev.K, dev.M))
    for j in np.nonzero(per_row)[0]:
        model = (
            radial(z[j]) if callable(radial) and not hasattr(radial, "cdf") else radial
        )
        pairs[j] = per_row[j] * model.annulus_fractions(faces)
    return pairs.ravel()


def ion_let_profile(particle, E0_MeV, mass_u=None):
    """(LET(z) keV/um, E_per_u(z) MeV/u) callables of the depth z (cm) for a proton or
    alpha of energy E0 entering SiC (NIST STAR, Bragg-additive SiC; CSDA residual
    energy, no straggling). Zero beyond the range."""
    from etna.core.stopping_power import linear_stopping, residual_energy_MeV

    A = {"proton": 1.00728, "alpha": 4.00151}[particle] if mass_u is None else mass_u

    def E_at(z):
        return np.asarray(residual_energy_MeV(particle, "sic", E0_MeV, z), dtype=float)

    def let(z):
        E = E_at(np.atleast_1d(z))
        out = np.zeros_like(E)
        ok = E > 1e-3
        out[ok] = linear_stopping(particle, "sic", E[ok]) * 1e-1  # MeV/cm -> keV/um
        return out

    def e_per_u(z):
        return E_at(np.atleast_1d(z)) / A

    return let, e_per_u


def whk_along(e_per_u, floor_MeV_u=1e-3):
    """Callable z -> WHK radial model at the local ion speed (``e_per_u`` from
    ion_let_profile)."""

    def model(z):
        return whk(beta_of(max(float(e_per_u(z)[0]), floor_MeV_u)))

    return model


def simulate_track(
    config=None,
    LET_keV_um=10.0,
    V_bias=50.0,
    radial=None,
    *,
    device=None,
    dark=None,
    mesh=None,
    region="epi",
    auger=True,
    pulse_sigma_s=None,
    rtol=1e-4,
    t_max=5e-6,
    thermal=None,
    **transient_kw,
):
    """One ion event: CCE, recombined fraction, ledger (``TrackResult``).

    ``config``: a DeviceConfig (None = the default stack) unless ``device`` (a
    TrackDevice) is given; ``dark`` = (x, xbc, system) from ``dark_state`` to reuse.
    ``thermal``: None (off), or a dict of ``LatticeHeat`` keywords (e.g.
    ``{"couple": True}``), or a LatticeHeat built on ``dark``'s system.
    """
    if device is None:
        from etna.api.device import DeviceConfig

        device = build_track_device(config or DeviceConfig(), mesh=mesh, auger=auger)
    if dark is None:
        dark = dark_state(device, V_bias)
    x_d, xbc, system = dark
    heat = None
    if thermal is not None:
        from .thermal import LatticeHeat

        heat = thermal if isinstance(thermal, LatticeHeat) else LatticeHeat(
            system, x_d, **dict(thermal)
        )  # fmt: skip
    pairs = track_pairs(device, LET_keV_um, radial, region)
    res = run_transient(system, x_d, xbc, pairs, pulse_sigma_s=pulse_sigma_s,
                        rtol=rtol, t_max=t_max, thermal=heat, **transient_kw)  # fmt: skip
    if heat is not None:
        res.extra["thermal"] = {
            "T_max_K": heat.T_max,
            "dT_max_K": heat.T_max - heat.T0,
            "heat_prompt_J": heat.heat_prompt,
            "heat_joule_J": heat.heat_joule,
            "heat_recomb_J": heat.heat_recomb,
            "heat_sink_J": heat.heat_sink,
            "stored_J": heat.stored(),
            "ledger_residual": heat.ledger_residual(),
        }
        if heat.couple:  # restore the isothermal coefficients for the next event
            system.set_temperature(None)
    return res


def cce_vs_let(
    config=None,
    LET_values=(0.1, 1.0, 10.0, 100.0, 1000.0),
    V_bias=50.0,
    radial=None,
    *,
    mesh=None,
    region="epi",
    auger=True,
    rtol=1e-4,
    **kw,
):
    """CCE / PHD table over ``LET_values`` (keV/um) on one device and dark state.

    Columns: LET_keV_um, CCE (pair-balance terminal charge / N_gen, TrackResult.cce),
    CCE_cathode (direct cathode flux sum, diagnostic), recombined_fraction, PHD (1 - CCE), peak_excess_cm3 (track-core mean
    density of the first cell at t = 0), ledger_residual, ledger_floor, converged,
    n_steps, t_end_ns, wall_s.
    """
    from etna.api.device import DeviceConfig

    dev = build_track_device(config or DeviceConfig(), mesh=mesh, auger=auger)
    dark = dark_state(dev, V_bias)
    rows = []
    for let in LET_values:
        t0 = time.perf_counter()
        r = simulate_track(device=dev, dark=dark, LET_keV_um=let, radial=radial,
                           region=region, rtol=rtol, V_bias=V_bias, **kw)  # fmt: skip
        rows.append(
            dict(
                LET_keV_um=float(let),
                CCE=r.cce,
                CCE_cathode=r.cce_cathode,
                recombined_fraction=r.recombined_fraction,
                PHD=1.0 - r.cce,
                peak_excess_cm3=r.peak_excess,
                ledger_residual=r.ledger_residual,
                ledger_floor=r.extra["ledger_floor"],
                converged=r.converged,
                n_steps=r.n_steps,
                t_end_ns=r.t_end * 1e9,
                wall_s=time.perf_counter() - t0,
            )
        )
    return pd.DataFrame(rows)


__all__ = [
    "MESH_DEFAULT",
    "MESH_FAST",
    "MESH_TINY",
    "TrackSystem",
    "build_track_device",
    "cce_vs_let",
    "ion_let_profile",
    "simulate_track",
    "track_pairs",
    "whk_along",
]

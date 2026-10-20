"""Monte Carlo coupling for 4H-SiC microdosimeter simulations.

Bridges external MC particle transport codes (Geant4/FLUKA) with the TCAD
charge collection pipeline.  Provides:

- CSV and ROOT event importers producing a standardized DataFrame
- Unit conversion (Geant4 mm/MeV -> project cm/keV)
- Per-event mesh charge profile mapping via ion_track_generation_2d
- Batch CCE(LET) lookup for fast ensemble processing (1000+ events/sec)
- Pulse height distribution histogram

Positions: both paths reduce an event to LET = sum(edep) / sv_thickness. The fast path
can weight the CCE by the event's lateral position (``process_mc_ensemble(...,
cce_lateral=...)``, audit-v7 B-11, 2026-10-20): the edep-weighted lateral centroid of
the event is mapped onto a 2D lateral CCE scan (charge_collection_2d.cce_lateral_scan,
normalised to its centre value). Without a map (the default) the CCE is
position-independent. The full path still places every track at ``x_ion_cm``. The default
2D device is laterally uniform (CCE(edge) / CCE(centre) = 1 by construction), so the map
only carries edge physics when it comes from a structure that has an edge (finite p+,
guard ring: alternative_structures).

Two processing paths:
    - **Full path** (events_to_charge_profiles): calls ion_track_generation_2d
      for each event, producing mesh-resolved 2D charge generation profiles.
      Suitable for detailed analysis of individual events.
    - **Fast path** (process_mc_ensemble): applies pre-computed CCE(LET) lookup
      table to each event's total deposited energy.  1000+ events in < 1 sec.

Units: Internally cm and keV.  Input units are configurable per loader.

Optional MC-truth columns (audit C-2): ``ekin_MeV`` -- kinetic energy (MeV, never
unit-converted) of the primary at the step, needed for the kinetic-energy kappa(E)
lookup in ``microdosimetry.tissue_equivalence_correction``; ``track_id`` -- if
present, the per-event kinetic energy is taken from the first row with
``track_id == 1`` (the Geant4 primary), else from the first row of the event.
``process_mc_ensemble`` returns it as ``event_kinetic_MeV`` (None without the column).

References:
    - Phase 21: single_particle.py (ion_track_generation_2d, load_cce_let_table)
    - Phase 22 research: MC event import architecture, CCE(LET) batch pattern
"""

import logging

import numpy as np
import pandas as pd

from etna.core.sic_material import SiC4H_Parameters

logger = logging.getLogger(__name__)

# E-h pair creation energy: single source of truth in SiC4H_Parameters (7.8 eV,
# measured 4H-SiC W-value). Audit M6: was hardcoded 8.4 eV here, ~7.7% high.
E_PAIR_SIC_EV = SiC4H_Parameters().E_pair_eV  # eV

# ---------------------------------------------------------------------------
# Unit conversion
# ---------------------------------------------------------------------------

UNIT_FACTORS = {
    "mm_to_cm": 0.1,
    "cm_to_cm": 1.0,
    "um_to_cm": 1e-4,
    "MeV_to_keV": 1e3,
    "keV_to_keV": 1.0,
    "eV_to_keV": 1e-3,
}


def convert_units(df, pos_unit="mm", energy_unit="MeV"):
    """Convert event DataFrame positions and energies to project units (cm, keV).

    Modifies x_cm, y_cm, z_cm and edep_keV columns in-place.

    Parameters
    ----------
    df : pd.DataFrame
        Event DataFrame with columns x_cm, y_cm, z_cm, edep_keV.
    pos_unit : str
        Source position unit: "mm", "cm", or "um".
    energy_unit : str
        Source energy unit: "MeV", "keV", or "eV".

    Returns
    -------
    df : pd.DataFrame
        Same DataFrame with converted values.

    Raises
    ------
    ValueError
        If pos_unit or energy_unit is not recognized.
    """
    pos_key = f"{pos_unit}_to_cm"
    energy_key = f"{energy_unit}_to_keV"

    if pos_key not in UNIT_FACTORS:
        raise ValueError(
            f"Unknown position unit '{pos_unit}'. " f"Supported: mm, cm, um"
        )
    if energy_key not in UNIT_FACTORS:
        raise ValueError(
            f"Unknown energy unit '{energy_unit}'. " f"Supported: MeV, keV, eV"
        )

    pos_factor = UNIT_FACTORS[pos_key]
    energy_factor = UNIT_FACTORS[energy_key]

    for col in ("x_cm", "y_cm", "z_cm"):
        if col in df.columns:
            df[col] = df[col] * pos_factor

    if "edep_keV" in df.columns:
        df["edep_keV"] = df["edep_keV"] * energy_factor

    return df


# ---------------------------------------------------------------------------
# CSV import
# ---------------------------------------------------------------------------


def load_mc_events_csv(filepath, column_map=None, pos_unit="cm", energy_unit="keV"):
    """Load Monte Carlo events from a CSV file.

    Parameters
    ----------
    filepath : str, pathlib.Path or file-like
        Path to the CSV file, or an open stream (anything ``pd.read_csv``
        accepts, e.g. ``io.BytesIO``).
    column_map : dict or None
        Maps standard column names to CSV column names.
        Keys are standard names (event_id, x_cm, y_cm, z_cm, edep_keV),
        values are the corresponding column names in the CSV file.
        Default: ``{"event_id": "event_id", "x_cm": "x", "y_cm": "y",
        "z_cm": "z", "edep_keV": "edep", "ekin_MeV": "ekin",
        "track_id": "track_id"}`` (the last two optional, see module docstring).
    pos_unit : str
        Position unit in the CSV file ("mm", "cm", "um").  Default: "cm".
    energy_unit : str
        Energy unit in the CSV file ("MeV", "keV", "eV").  Default: "keV".

    Returns
    -------
    events : pd.DataFrame
        Standardized DataFrame with columns: event_id, x_cm, y_cm, z_cm,
        edep_keV.
    """
    default_map = {
        "event_id": "event_id",
        "x_cm": "x",
        "y_cm": "y",
        "z_cm": "z",
        "edep_keV": "edep",
        "ekin_MeV": "ekin",
        "track_id": "track_id",
    }
    col_map = default_map if column_map is None else column_map

    df = pd.read_csv(filepath)

    # Rename: CSV column name -> standard name
    rename = {
        csv_col: std_name
        for std_name, csv_col in col_map.items()
        if csv_col in df.columns
    }
    df = df.rename(columns=rename)

    # Apply unit conversion if needed
    if pos_unit != "cm" or energy_unit != "keV":
        convert_units(df, pos_unit=pos_unit, energy_unit=energy_unit)

    n_events = df["event_id"].nunique()
    n_rows = len(df)
    logger.info("Loaded %d rows (%d unique events) from %s", n_rows, n_events, filepath)

    return df


# ---------------------------------------------------------------------------
# ROOT import
# ---------------------------------------------------------------------------


def list_root_trees(filepath):
    """List TTree names in a ROOT file.

    Kept as public API (no internal caller): the discovery step for
    ``load_mc_events_root(tree_name=...)`` when a Geant4 file does not use "Hits".

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to ROOT file.

    Returns
    -------
    tree_names : list of str
        Names of TTree objects in the file (without cycle numbers).

    Raises
    ------
    FileNotFoundError
        If the ROOT file does not exist.
    """
    try:
        import uproot
    except ImportError:
        raise ImportError(
            "uproot is required for ROOT file support. "
            "Install with: uv pip install uproot>=5.6"
        )

    try:
        with uproot.open(filepath) as f:
            all_keys = f.keys()
            # Filter to TTree objects only
            tree_names = []
            for key in all_keys:
                # Strip cycle number (e.g., "Hits;1" -> "Hits")
                name = key.split(";")[0]
                obj = f[key]
                # Check if it's a TTree (has .keys() method for branches)
                if hasattr(obj, "num_entries"):
                    if name not in tree_names:
                        tree_names.append(name)
            return tree_names
    except FileNotFoundError:
        raise FileNotFoundError(
            f"ROOT file not found: {filepath}. "
            "Check the file path and ensure the file exists."
        )


def load_mc_events_root(
    filepath,
    tree_name="Hits",
    branch_map=None,
    pos_unit="mm",
    energy_unit="MeV",
):
    """Load Monte Carlo events from a Geant4 ROOT file via uproot.

    Parameters
    ----------
    filepath : str or pathlib.Path
        Path to ROOT file.
    tree_name : str
        Name of TTree in ROOT file.  Default: "Hits".
    branch_map : dict or None
        Maps standard column names to ROOT branch names.
        Keys are standard names (event_id, x_cm, y_cm, z_cm, edep_keV),
        values are ROOT branch names.
        Default: ``{"event_id": "EventID", "x_cm": "PosX", "y_cm": "PosY",
        "z_cm": "PosZ", "edep_keV": "Edep"}``.
    pos_unit : str
        Position unit in ROOT file ("mm", "cm", "um").  Default: "mm"
        (Geant4 default).
    energy_unit : str
        Energy unit in ROOT file ("MeV", "keV", "eV").  Default: "MeV"
        (Geant4 default).

    Returns
    -------
    events : pd.DataFrame
        Standardized DataFrame with columns: event_id, x_cm, y_cm, z_cm,
        edep_keV.  Optional: particle, ekin_MeV, track_id (if in branch_map and
        available; ekin_MeV must already be in MeV).
    """
    try:
        import uproot
    except ImportError:
        raise ImportError(
            "uproot is required for ROOT file support. "
            "Install with: uv pip install uproot>=5.6"
        )

    default_map = {
        "event_id": "EventID",
        "x_cm": "PosX",
        "y_cm": "PosY",
        "z_cm": "PosZ",
        "edep_keV": "Edep",
    }
    bmap = default_map if branch_map is None else branch_map

    with uproot.open(filepath) as f:
        tree = f[tree_name]
        available_branches = set(tree.keys())

        # Determine which branches to read
        branches_to_read = {}
        for std_name, root_name in bmap.items():
            if root_name in available_branches:
                branches_to_read[std_name] = root_name
            elif std_name in ("particle", "ekin_MeV", "track_id"):
                # Optional branch -- skip if not present
                logger.debug("Optional branch '%s' not found in tree", root_name)
            else:
                raise KeyError(
                    f"Required branch '{root_name}' (for '{std_name}') "
                    f"not found in tree '{tree_name}'. "
                    f"Available branches: {sorted(available_branches)}"
                )

        arrays = tree.arrays(list(branches_to_read.values()), library="np")

    df = pd.DataFrame(
        {
            std_name: arrays[root_name]
            for std_name, root_name in branches_to_read.items()
        }
    )

    # Apply unit conversion (Geant4 default: mm, MeV)
    if pos_unit != "cm" or energy_unit != "keV":
        convert_units(df, pos_unit=pos_unit, energy_unit=energy_unit)

    n_events = df["event_id"].nunique()
    n_rows = len(df)
    logger.info(
        "Loaded %d rows (%d unique events) from %s:%s",
        n_rows,
        n_events,
        filepath,
        tree_name,
    )

    return df


# ---------------------------------------------------------------------------
# Mesh charge profile mapping (full path)
# ---------------------------------------------------------------------------


def events_to_charge_profiles(
    events_df,
    device_info,
    sv_thickness_um=10.0,
    x_ion_cm=0.0,
    track_sigma_cm=1e-4,
    E_pair_eV=E_PAIR_SIC_EV,
):
    """Convert MC events to 2D mesh charge generation profiles.

    Public pipeline step (no internal caller; listed in the module docstring): Geant4
    events -> per-event 2D generation for the single-particle solver.

    This is the "full path" for mesh-resolved profiles.  Each event's total
    deposited energy is converted to an effective LET, then
    ``ion_track_generation_2d`` produces a spatially-resolved charge generation
    profile on the 2D devsim mesh.

    For fast batch processing of 1000+ events (without per-event mesh
    generation), use :func:`process_mc_ensemble` instead.

    Parameters
    ----------
    events_df : pd.DataFrame
        Standardized event DataFrame (step-level or event-level).
    device_info : dict
        2D device info from create_2d_dd_device.
    sv_thickness_um : float
        Sensitive volume thickness for LET computation (um).
    x_ion_cm : float
        Lateral position of ion track (cm).
    track_sigma_cm : float
        Gaussian radial width of ion track (cm).
    E_pair_eV : float
        Electron-hole pair creation energy (eV).

    Returns
    -------
    profiles : list of dict
        One dict per event with keys: event_id (int), LET_keV_um (float),
        generation (ndarray), Q_generated_C_per_cm (float).
    """
    from etna.core.single_particle import ion_track_generation_2d

    # Sum energy per event
    event_totals = events_df.groupby("event_id")["edep_keV"].sum()

    profiles = []
    n_events = len(event_totals)
    for i, (event_id, total_edep) in enumerate(event_totals.items()):
        LET_keV_um = total_edep / sv_thickness_um

        generation, Q_gen = ion_track_generation_2d(
            device_info,
            LET_keV_um,
            x_ion_cm=x_ion_cm,
            track_sigma_cm=track_sigma_cm,
            E_pair_eV=E_pair_eV,
        )

        profiles.append(
            {
                "event_id": int(event_id),
                "LET_keV_um": float(LET_keV_um),
                "generation": generation,
                "Q_generated_C_per_cm": float(Q_gen),
            }
        )

        if (i + 1) % 100 == 0:
            logger.info("Processed %d / %d events for charge profiles", i + 1, n_events)

    logger.info("Generated charge profiles for %d events", n_events)
    return profiles


# ---------------------------------------------------------------------------
# Fast batch CCE lookup (fast path)
# ---------------------------------------------------------------------------


def lateral_cce_factor(r_cm, cce_lateral, outside="edge"):
    """Relative lateral collection f(r) = CCE(r) / CCE(0) of a 2D lateral scan.

    ``cce_lateral``: a :func:`etna.core.charge_collection_2d.cce_lateral_scan` result
    (or any dict with ``x_positions_cm`` from the SV centre outwards and
    ``cce_values``); ``r_cm`` (array) distances from the SV centre, linearly
    interpolated. Beyond the last scan position: ``outside="edge"`` uses the last value
    (clamped), ``"zero"`` gives 0 (no collection outside the scanned device), ``"raise"``
    raises ValueError. Returns (f, n_outside).
    """
    x = np.asarray(cce_lateral["x_positions_cm"], dtype=float)
    c = np.asarray(cce_lateral["cce_values"], dtype=float)
    if x.ndim != 1 or x.size < 2 or x.shape != c.shape:
        raise ValueError(
            "cce_lateral needs 1D x_positions_cm / cce_values of equal length"
        )
    if np.any(np.diff(x) <= 0) or x[0] != 0.0:
        raise ValueError("cce_lateral x_positions_cm must start at 0 and increase")
    if not (np.all(np.isfinite(c)) and c[0] > 0):
        raise ValueError("cce_lateral cce_values must be finite with CCE(0) > 0")
    if outside not in ("edge", "zero", "raise"):
        raise ValueError(f"outside must be 'edge', 'zero' or 'raise' (got {outside!r})")
    r = np.abs(np.asarray(r_cm, dtype=float))
    beyond = r > x[-1]
    n_out = int(np.sum(beyond))
    if n_out and outside == "raise":
        raise ValueError(
            f"{n_out} events lie beyond the lateral CCE scan ({x[-1] * 1e4:g} um)"
        )
    f = np.interp(r, x, c / c[0])
    if outside == "zero":
        f = np.where(beyond, 0.0, f)
    return f, n_out


def event_lateral_distance(events_df, sv_center_cm=(0.0, 0.0), metric="chebyshev"):
    """Per-event lateral distance (cm) of the edep-weighted hit centroid from the SV
    centre, indexed by event_id (sorted, as ``groupby``).

    Lateral plane = (x_cm, y_cm); depth = z_cm (the loaders' convention). ``metric``:
    ``"chebyshev"`` max(|dx|, |dy|) (square SV: distance towards the nearest side, so a
    1D centre-to-edge scan applies along both axes), ``"radial"`` hypot(dx, dy) (round
    SV), ``"x"`` |dx| (strip SV, edges along x only). Events with zero total edep get
    the plain mean position.
    """
    if metric not in ("chebyshev", "radial", "x"):
        raise ValueError(
            f"metric must be 'chebyshev', 'radial' or 'x' (got {metric!r})"
        )
    cols = ["x_cm", "y_cm"] if metric != "x" else ["x_cm"]
    missing = [c for c in cols if c not in events_df.columns]
    if missing:
        raise ValueError(f"lateral CCE weighting needs position columns {missing}")
    w = events_df["edep_keV"].clip(lower=0.0)
    g = events_df.assign(_w=w)
    out = {}
    tot = g.groupby("event_id")["_w"].sum()
    for i, c in enumerate(cols):
        num = (g[c] * g["_w"]).groupby(g["event_id"]).sum()
        mean = g.groupby("event_id")[c].mean()
        has_w = tot.values > 0
        safe = np.where(has_w, tot.values, 1.0)
        centroid = np.where(has_w, num.values / safe, mean.values)
        out[c] = centroid - float(sv_center_cm[i])
    if metric == "x":
        d = np.abs(out["x_cm"])
    elif metric == "radial":
        d = np.hypot(out["x_cm"], out["y_cm"])
    else:
        d = np.maximum(np.abs(out["x_cm"]), np.abs(out["y_cm"]))
    return pd.Series(d, index=tot.index)


def process_mc_ensemble(
    events_df,
    cce_interp,
    sv_thickness_um=10.0,
    cce_lateral=None,
    sv_center_cm=(0.0, 0.0),
    lateral_metric="chebyshev",
    lateral_outside="edge",
):
    """Process MC events through CCE(LET) lookup table.

    Groups step-level data by event_id, computes effective LET per event,
    applies the pre-computed CCE interpolation function, and returns
    collected energy per event.

    Lateral (edge) weighting, optional (audit-v7 B-11): with ``cce_lateral`` (a
    :func:`etna.core.charge_collection_2d.cce_lateral_scan` result) each event's CCE is
    CCE(LET) x f(r), f = CCE_2D(r) / CCE_2D(0) (lateral_cce_factor), r the distance of
    its edep-weighted lateral centroid from ``sv_center_cm`` (event_lateral_distance,
    ``lateral_metric``). Normalising to the centre keeps the CCE(LET) table (a centre /
    1D quantity) from being counted twice. Assumption: separable CCE(LET, r) = CCE(LET)
    f(r), with f from the scan's own (alpha-stripe) generation; the centroid stands for
    the whole track (tracks crossing the edge are averaged, not split).

    Parameters
    ----------
    events_df : pd.DataFrame
        Standardized event DataFrame (step-level or event-level).
    cce_interp : callable
        CCE interpolation function from ``load_cce_let_table()``.
        Signature: cce_interp(LET_keV_um) -> CCE. If it carries ``let_min`` /
        ``let_max`` (load_cce_let_table sets them), events outside that range are
        counted in ``n_clamped`` and logged: their CCE is the table edge value.
    sv_thickness_um : float
        Sensitive volume thickness for LET computation (um).
    cce_lateral : dict, optional
        Lateral CCE scan (``x_positions_cm``, ``cce_values``); None = no lateral
        weighting (position-independent CCE, the behaviour before 2026-10-20).
    sv_center_cm : (float, float)
        SV centre in the event (x_cm, y_cm) frame. Default (0, 0).
    lateral_metric : {"chebyshev", "radial", "x"}
        Lateral distance of event_lateral_distance. Default "chebyshev" (square SV).
    lateral_outside : {"edge", "zero", "raise"}
        Events beyond the last scan position (lateral_cce_factor). Default "edge"
        (clamped, counted in ``n_outside_lateral`` and logged).

    Returns
    -------
    result : dict
        - event_energies_keV : ndarray -- total deposited energy per event
        - event_LET_keV_um : ndarray -- effective LET per event
        - event_CCE : ndarray -- charge collection efficiency per event
        - event_collected_keV : ndarray -- collected energy per event
        - n_events : int -- number of events processed
        - n_zero_energy : int -- number of events with zero/negative energy
        - n_clamped : int -- events with LET outside the CCE table (edge CCE used);
          0 when cce_interp carries no range
        - event_lateral_factor : ndarray -- f(r) per event (1 without cce_lateral)
        - event_lateral_distance_cm : ndarray or None -- r per event
        - n_outside_lateral : int -- events beyond the lateral scan (0 without it)
        - event_kinetic_MeV : ndarray or None -- primary kinetic energy per event
          (same order/mask as event_energies_keV) if the events carry ``ekin_MeV``
          (audit C-2: the variable for the kappa(E) lookup), else None
    """
    # Sum energy per event
    event_totals = events_df.groupby("event_id")["edep_keV"].sum()
    energies = event_totals.values.astype(float)
    ekin = _event_kinetic_energy(events_df, event_totals.index)

    # Filter zero-energy events
    mask_positive = energies > 0
    n_zero = int(np.sum(~mask_positive))
    if n_zero > 0:
        logger.info(
            "Filtered %d zero/negative-energy events (%.1f%% of total)",
            n_zero,
            100.0 * n_zero / len(energies),
        )

    energies_pos = energies[mask_positive]

    # Compute effective LET
    LET = energies_pos / sv_thickness_um  # keV/um

    # Warn if LET outside typical range
    if len(LET) > 0:
        if np.any(LET < 0.01) or np.any(LET > 10000):
            n_low = int(np.sum(LET < 0.01))
            n_high = int(np.sum(LET > 10000))
            logger.warning(
                "LET outside typical range [0.01, 10000] keV/um: " "%d below, %d above",
                n_low,
                n_high,
            )

    # CCE table range: the lookup clamps outside it (np.interp edge value)
    n_clamped = 0
    lo = getattr(cce_interp, "let_min", None)
    hi = getattr(cce_interp, "let_max", None)
    if lo is not None and hi is not None and len(LET) > 0:
        n_clamped = int(np.sum((LET < lo) | (LET > hi)))
        if n_clamped:
            logger.warning(
                "%d of %d events have LET outside the CCE table [%g, %g] keV/um: "
                "their CCE is the table edge value (clamped)",
                n_clamped,
                len(LET),
                lo,
                hi,
            )

    # Vectorized CCE lookup
    vectorized = getattr(cce_interp, "vectorized", False)  # load_cce_let_table's
    CCE = (
        np.asarray(cce_interp(LET), dtype=float)
        if vectorized
        else np.array([cce_interp(let) for let in LET])
    )

    # Lateral (edge) weighting of the CCE by the event position (audit-v7 B-11)
    f_lat = np.ones(len(LET))
    r_evt = None
    n_out = 0
    if cce_lateral is not None:
        r_all = event_lateral_distance(events_df, sv_center_cm, lateral_metric)
        r_evt = r_all.reindex(event_totals.index).to_numpy()[mask_positive]
        f_lat, n_out = lateral_cce_factor(r_evt, cce_lateral, outside=lateral_outside)
        if n_out:
            logger.warning(
                "%d of %d events lie beyond the lateral CCE scan (%g um): %s",
                n_out,
                len(LET),
                float(np.max(cce_lateral["x_positions_cm"])) * 1e4,
                "edge value used" if lateral_outside == "edge" else "CCE set to 0",
            )
        CCE = CCE * f_lat

    # Collected energy
    E_collected = CCE * energies_pos

    return {
        "event_energies_keV": energies_pos,
        "event_LET_keV_um": LET,
        "event_CCE": CCE,
        "event_collected_keV": E_collected,
        "n_events": int(len(energies_pos)),
        "n_zero_energy": n_zero,
        "n_clamped": n_clamped,
        "event_lateral_factor": f_lat,
        "event_lateral_distance_cm": r_evt,
        "n_outside_lateral": n_out,
        "event_kinetic_MeV": None if ekin is None else ekin[mask_positive],
    }


def _event_kinetic_energy(events_df, event_index):
    """Per-event primary kinetic energy (MeV) aligned with ``event_index``, or None.

    First row with ``track_id == 1`` if a track_id column exists (Geant4 primary),
    else the first row of the event (file order). NaN for events without such a row.
    """
    if "ekin_MeV" not in events_df.columns:
        return None
    df = events_df
    if "track_id" in df.columns:
        df = df[df["track_id"] == 1]
    first = df.groupby("event_id", sort=True)["ekin_MeV"].first()
    return first.reindex(event_index).to_numpy(dtype=float)


# ---------------------------------------------------------------------------
# Pulse height distribution
# ---------------------------------------------------------------------------


def pulse_height_distribution(
    collected_energies_keV, n_bins=200, e_min=None, e_max=None
):
    """Build pulse height distribution from collected energies.

    Parameters
    ----------
    collected_energies_keV : array_like
        Collected energy per event (keV).
    n_bins : int
        Number of histogram bins.  Default: 200.
    e_min : float or None
        Minimum energy for histogram range (keV).  Auto-detected if None
        (0.5 * min of positive energies).
    e_max : float or None
        Maximum energy for histogram range (keV).  Auto-detected if None
        (1.5 * max of positive energies).

    Returns
    -------
    result : dict
        - bin_centers_keV : ndarray -- geometric mean of bin edges (keV)
        - counts : ndarray -- event counts per bin
        - bin_edges_keV : ndarray -- bin edge energies (keV)
    """
    E = np.asarray(collected_energies_keV, dtype=float)
    E = E[E > 0]

    if len(E) == 0:
        logger.warning("No positive-energy events for pulse height distribution")
        edges = np.logspace(0, 3, n_bins + 1)
        return {
            "bin_centers_keV": np.sqrt(edges[:-1] * edges[1:]),
            "counts": np.zeros(n_bins, dtype=int),
            "bin_edges_keV": edges,
        }

    if e_min is None:
        e_min = E.min() * 0.5
    if e_max is None:
        e_max = E.max() * 1.5

    bin_edges = np.logspace(np.log10(e_min), np.log10(e_max), n_bins + 1)
    counts, _ = np.histogram(E, bins=bin_edges)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # geometric mean

    return {
        "bin_centers_keV": bin_centers,
        "counts": counts,
        "bin_edges_keV": bin_edges,
    }

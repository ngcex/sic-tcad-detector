"""Monte Carlo coupling for 4H-SiC microdosimeter simulations.

Bridges external MC particle transport codes (Geant4/FLUKA) with the TCAD
charge collection pipeline.  Provides:

- CSV and ROOT event importers producing a standardized DataFrame
- Unit conversion (Geant4 mm/MeV -> project cm/keV)
- Per-event mesh charge profile mapping via ion_track_generation_2d
- Batch CCE(LET) lookup for fast ensemble processing (1000+ events/sec)
- Pulse height distribution histogram

Two processing paths:
    - **Full path** (events_to_charge_profiles): calls ion_track_generation_2d
      for each event, producing mesh-resolved 2D charge generation profiles.
      Suitable for detailed analysis of individual events.
    - **Fast path** (process_mc_ensemble): applies pre-computed CCE(LET) lookup
      table to each event's total deposited energy.  1000+ events in < 1 sec.

Units: Internally cm and keV.  Input units are configurable per loader.

References:
    - Phase 21: single_particle.py (ion_track_generation_2d, load_cce_let_table)
    - Phase 22 research: MC event import architecture, CCE(LET) batch pattern
"""

import logging

import numpy as np
import pandas as pd

logger = logging.getLogger(__name__)

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
    filepath : str or pathlib.Path
        Path to the CSV file.
    column_map : dict or None
        Maps standard column names to CSV column names.
        Keys are standard names (event_id, x_cm, y_cm, z_cm, edep_keV),
        values are the corresponding column names in the CSV file.
        Default: ``{"event_id": "event_id", "x_cm": "x", "y_cm": "y",
        "z_cm": "z", "edep_keV": "edep"}``.
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
        edep_keV.  Optional: particle (if in branch_map and available).
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
            elif std_name == "particle":
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
    E_pair_eV=8.4,
):
    """Convert MC events to 2D mesh charge generation profiles.

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
    from src.single_particle import ion_track_generation_2d

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


def process_mc_ensemble(events_df, cce_interp, sv_thickness_um=10.0):
    """Process MC events through CCE(LET) lookup table.

    Groups step-level data by event_id, computes effective LET per event,
    applies the pre-computed CCE interpolation function, and returns
    collected energy per event.

    Parameters
    ----------
    events_df : pd.DataFrame
        Standardized event DataFrame (step-level or event-level).
    cce_interp : callable
        CCE interpolation function from ``load_cce_let_table()``.
        Signature: cce_interp(LET_keV_um) -> CCE.
    sv_thickness_um : float
        Sensitive volume thickness for LET computation (um).

    Returns
    -------
    result : dict
        - event_energies_keV : ndarray -- total deposited energy per event
        - event_LET_keV_um : ndarray -- effective LET per event
        - event_CCE : ndarray -- charge collection efficiency per event
        - event_collected_keV : ndarray -- collected energy per event
        - n_events : int -- number of events processed
        - n_zero_energy : int -- number of events with zero/negative energy
    """
    # Sum energy per event
    event_totals = events_df.groupby("event_id")["edep_keV"].sum()
    energies = event_totals.values.astype(float)

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

    # Vectorized CCE lookup
    CCE = np.array([cce_interp(let) for let in LET])

    # Collected energy
    E_collected = CCE * energies_pos

    return {
        "event_energies_keV": energies_pos,
        "event_LET_keV_um": LET,
        "event_CCE": CCE,
        "event_collected_keV": E_collected,
        "n_events": int(len(energies_pos)),
        "n_zero_energy": n_zero,
    }


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

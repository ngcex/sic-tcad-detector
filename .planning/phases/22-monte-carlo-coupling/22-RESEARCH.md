# Phase 22: Monte Carlo Coupling - Research

**Researched:** 2026-03-31
**Domain:** MC energy deposition import (CSV + ROOT), event-to-mesh mapping, CCE(LET) lookup application, pulse height distribution
**Confidence:** HIGH

## Summary

Phase 22 bridges external Monte Carlo particle transport codes (Geant4/FLUKA) with the TCAD charge collection pipeline built in Phases 19-21. The core workflow is: (1) import per-event energy deposition data from CSV or ROOT files, (2) convert each event's energy deposition into a charge generation profile on the 2D devsim mesh, and (3) apply the pre-computed CCE(LET) lookup table to map deposited energy to collected charge for thousands of events, producing a pulse height distribution.

The existing `single_particle.py` module (Phase 21) provides `load_cce_let_table()` which returns a `cce_interp(LET_keV_um) -> CCE` function using log-linear interpolation, and `ion_track_generation_2d()` for converting LET to mesh-resolved charge profiles. The CCE(LET) table is stored as JSON with geometry metadata. The key new capability is: (1) parsers for CSV and Geant4 ROOT event files, (2) per-step or per-event energy-to-LET conversion, and (3) a batch processing loop that applies CCE interpolation to produce collected charge per event, yielding a pulse height histogram.

For the ROOT file import (MCCP-02), `uproot` (>=5.6, already declared as a v3.0 dependency) reads TTree branches as numpy arrays without requiring C++ ROOT. The Geant4 TTree naming conventions from the INFN-LNS group are unknown (flagged as a blocker in STATE.md), so the ROOT importer must accept configurable branch names. For CSV import (MCCP-01), a simple pandas `read_csv` with column name mapping suffices. The mesh mapping (MCCP-03) reuses `ion_track_generation_2d()` for vertical tracks, or for events with step-level spatial data, sums per-step contributions. The batch CCE application (MCCP-04) is a vectorized numpy operation: for each event, compute total energy deposited in the SV, convert to LET, look up CCE, and compute collected charge = CCE _ deposited_energy / E_pair _ q.

**Primary recommendation:** Create a `mc_coupling.py` module with: (1) `load_mc_events_csv()` and `load_mc_events_root()` parsers returning a standardized event DataFrame, (2) `events_to_charge_profiles()` for mesh mapping, (3) `process_mc_ensemble()` for batch CCE lookup producing a pulse height distribution. Keep ROOT branch names configurable with sensible defaults.

<phase_requirements>

## Phase Requirements

| ID      | Description                                                                                        | Research Support                                                                                                                                                                                                                                                       |
| ------- | -------------------------------------------------------------------------------------------------- | ---------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| MCCP-01 | Import energy deposition from CSV (position + energy columns) for any ion species                  | pandas `read_csv` with configurable column mapping. Standardize to DataFrame with columns: event_id, x_cm, y_cm, z_cm, edep_keV. Column name aliases via dict parameter.                                                                                               |
| MCCP-02 | Import energy deposition from Geant4 ROOT files using uproot (TTree with position/energy branches) | `uproot.open(path)[tree_name].arrays(branch_list, library="np")`. Branch names configurable (blocker: INFN-LNS naming unknown). Default branch names from Geant4 analysis manager conventions.                                                                         |
| MCCP-03 | Convert MC energy deposition events to charge generation profiles on 2D devsim mesh                | Two modes: (a) fast path -- sum edep per event to get total LET, apply existing `ion_track_generation_2d()` for center-track profile; (b) full path -- per-step spatial mapping onto mesh nodes for resolved tracks. Fast path sufficient for CCE(LET) table approach. |
| MCCP-04 | Process 1000+ MC events through CCE(LET) lookup to build pulse height distribution                 | Vectorized: LET_per_event = total_edep_keV / thickness_um. CCE = cce_interp(LET). E_collected = CCE \* total_edep. Histogram with np.histogram on log-spaced bins. Must handle 1000+ events in <1 second.                                                              |

</phase_requirements>

## Standard Stack

### Core

| Library    | Version | Purpose                                          | Why Standard                                                      |
| ---------- | ------- | ------------------------------------------------ | ----------------------------------------------------------------- |
| uproot     | >=5.6   | Read Geant4 ROOT TTree files as numpy arrays     | Already declared v3.0 dependency; pure Python, no C++ ROOT needed |
| numpy      | >=1.24  | Array ops, histogram, vectorized CCE lookup      | Already in stack                                                  |
| pandas     | >=2.0   | CSV parsing, event DataFrame, tabular operations | Already in stack                                                  |
| matplotlib | >=3.7   | Pulse height distribution plots                  | Already in stack                                                  |

### Supporting

| Library | Version  | Purpose                                 | When to Use                         |
| ------- | -------- | --------------------------------------- | ----------------------------------- |
| json    | (stdlib) | Load CCE(LET) table from Phase 21       | Every run -- `load_cce_let_table()` |
| pathlib | (stdlib) | File path handling                      | File I/O                            |
| logging | (stdlib) | Progress reporting for batch processing | Always                              |

### Alternatives Considered

| Instead of               | Could Use                | Tradeoff                                                                                                                                                                       |
| ------------------------ | ------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------ |
| uproot for ROOT          | PyROOT (C++ binding)     | PyROOT requires C++ ROOT installation; uproot is pure Python. Use uproot.                                                                                                      |
| pandas for CSV           | numpy.loadtxt            | pandas handles headers, mixed types, missing values better. Use pandas.                                                                                                        |
| Per-event TCAD transient | CCE(LET) lookup table    | Full transient per event would take ~30 sec \* 1000 events = 8+ hours. Lookup table gives result in <1 sec. Use lookup table.                                                  |
| Per-step spatial mapping | Event-total LET approach | Step-level mapping preserves track structure but is unnecessary when using CCE(LET) table (which assumes vertical center track). Implement both; default to fast LET approach. |

**Installation:**

```bash
uv pip install uproot>=5.6
```

Note: uproot is already declared as a v3.0 dependency. Verify it is installed.

## Architecture Patterns

### Recommended Project Structure

```
src/
  mc_coupling.py          # NEW: MC event import, mesh mapping, batch CCE
  single_particle.py      # Phase 21: CCE(LET) table, load_cce_let_table()
  charge_collection_2d.py # Phase 20: create_2d_dd_device, integrate_over_mesh_2d
  device2d.py             # Phase 19: 2D mesh generation
notebooks/
  17_mc_coupling.ipynb    # NEW: publication notebook for Phase 22
```

### Pattern 1: Standardized Event DataFrame

**What:** All importers (CSV, ROOT) produce a common DataFrame format regardless of source.
**When to use:** Always -- decouple parsing from processing.
**Example:**

```python
# Standardized event format after import
# Columns: event_id (int), x_cm (float), y_cm (float), z_cm (float),
#          edep_keV (float), particle (str, optional)
#
# Each row is one energy deposition step.
# Multiple rows per event_id = step-level data.
# One row per event_id = event-total data.

def load_mc_events_csv(filepath, column_map=None):
    """Load MC events from CSV with configurable column mapping.

    Parameters
    ----------
    filepath : str or Path
        Path to CSV file.
    column_map : dict or None
        Maps CSV column names to standard names.
        Default: {"event_id": "event_id", "x": "x_cm", "y": "y_cm",
                  "z": "z_cm", "edep": "edep_keV"}

    Returns
    -------
    events : pd.DataFrame
        Standardized event DataFrame.
    """
    default_map = {
        "event_id": "event_id",
        "x": "x_cm",
        "y": "y_cm",
        "z": "z_cm",
        "edep": "edep_keV",
    }
    col_map = default_map if column_map is None else column_map
    df = pd.read_csv(filepath)
    df = df.rename(columns={v: k for k, v in col_map.items()
                            if v in df.columns})
    return df
```

### Pattern 2: Configurable ROOT Branch Names

**What:** ROOT importer accepts branch name mapping because Geant4 TTree naming varies by group/simulation.
**When to use:** Always for ROOT import -- the INFN-LNS naming is unknown.
**Example:**

```python
import uproot
import numpy as np

def load_mc_events_root(filepath, tree_name="Hits", branch_map=None):
    """Load MC events from Geant4 ROOT file via uproot.

    Parameters
    ----------
    filepath : str or Path
        Path to ROOT file.
    tree_name : str
        Name of TTree in ROOT file. Default: "Hits".
    branch_map : dict or None
        Maps standard names to ROOT branch names.
        Default: {"event_id": "EventID", "x_cm": "PosX",
                  "y_cm": "PosY", "z_cm": "PosZ",
                  "edep_keV": "Edep"}

    Returns
    -------
    events : pd.DataFrame
        Standardized event DataFrame.
    """
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
        branches = list(bmap.values())
        arrays = tree.arrays(branches, library="np")

    df = pd.DataFrame({
        std_name: arrays[root_name]
        for std_name, root_name in bmap.items()
    })
    return df
```

### Pattern 3: Fast Vectorized CCE Application

**What:** Process entire event ensemble without per-event TCAD simulation.
**When to use:** For MCCP-04 batch processing.
**Example:**

```python
def process_mc_ensemble(events_df, cce_interp, sv_thickness_um=10.0,
                        E_pair_eV=8.4):
    """Process MC events through CCE(LET) lookup.

    Parameters
    ----------
    events_df : pd.DataFrame
        Standardized event DataFrame (step-level or event-level).
    cce_interp : callable
        CCE interpolation function from load_cce_let_table().
    sv_thickness_um : float
        SV thickness for LET calculation (um).
    E_pair_eV : float
        Electron-hole pair creation energy (eV).

    Returns
    -------
    result : dict
        "event_energies_keV": total deposited energy per event (array)
        "event_LET_keV_um": effective LET per event (array)
        "event_CCE": CCE per event (array)
        "event_collected_keV": collected energy per event (array)
    """
    # Sum energy per event
    event_totals = events_df.groupby("event_id")["edep_keV"].sum()

    # Effective LET = total_edep / thickness
    LET = event_totals.values / sv_thickness_um  # keV/um

    # Vectorized CCE lookup
    CCE = np.array([cce_interp(let) for let in LET])

    # Collected energy
    E_collected = CCE * event_totals.values

    return {
        "event_energies_keV": event_totals.values,
        "event_LET_keV_um": LET,
        "event_CCE": CCE,
        "event_collected_keV": E_collected,
    }
```

### Pattern 4: Pulse Height Distribution

**What:** Histogram of collected energies across all events.
**When to use:** Final output of MCCP-04.
**Example:**

```python
def pulse_height_distribution(collected_energies_keV, n_bins=200,
                               e_min=None, e_max=None):
    """Build pulse height distribution from collected energies.

    Parameters
    ----------
    collected_energies_keV : array_like
        Collected energy per event (keV).
    n_bins : int
        Number of histogram bins.
    e_min, e_max : float or None
        Energy range (keV). Auto-detected if None.

    Returns
    -------
    bin_centers : ndarray
        Bin center energies (keV).
    counts : ndarray
        Event counts per bin.
    bin_edges : ndarray
        Bin edge energies (keV).
    """
    E = np.asarray(collected_energies_keV)
    E = E[E > 0]  # exclude zero-energy events

    if e_min is None:
        e_min = E.min() * 0.5
    if e_max is None:
        e_max = E.max() * 1.5

    bin_edges = np.logspace(np.log10(e_min), np.log10(e_max), n_bins + 1)
    counts, _ = np.histogram(E, bins=bin_edges)
    bin_centers = np.sqrt(bin_edges[:-1] * bin_edges[1:])  # geometric mean

    return bin_centers, counts, bin_edges
```

### Anti-Patterns to Avoid

- **Running TCAD transient per MC event:** The CCE(LET) lookup table exists precisely to avoid this. Each transient takes ~30 sec; 1000 events would take 8+ hours. Use the lookup table.
- **Hardcoding ROOT branch names:** Geant4 output naming varies by user/group. Always accept configurable branch/tree names.
- **Ignoring event boundaries:** MC output files may have step-level data (multiple rows per event). Always group by event_id before computing per-event quantities.
- **Assuming units:** Geant4 default units are mm and MeV. The project uses cm and keV. Unit conversion must be explicit and validated.

## Don't Hand-Roll

| Problem                  | Don't Build                 | Use Instead                                    | Why                                                                                    |
| ------------------------ | --------------------------- | ---------------------------------------------- | -------------------------------------------------------------------------------------- |
| ROOT file reading        | Custom binary parser        | uproot                                         | ROOT format is complex; uproot handles compression, TTree serialization, jagged arrays |
| CSV parsing with headers | Manual string splitting     | pandas.read_csv                                | Handles quoting, missing values, type inference, encoding                              |
| CCE interpolation        | New interpolation code      | `load_cce_let_table()` from single_particle.py | Already validated, handles NaN filtering, log-linear interp                            |
| Histogram binning        | Manual bin counting         | numpy.histogram                                | Handles edge cases, overflow, underflow correctly                                      |
| Unit conversion          | Ad-hoc per-function factors | Centralized conversion constants               | Geant4 uses mm/MeV, project uses cm/keV; one place for factors                         |

**Key insight:** This phase is primarily data plumbing -- importing, converting, and applying an existing lookup table. The physics complexity was solved in Phase 21. The challenge here is robust I/O with configurable formats.

## Common Pitfalls

### Pitfall 1: Geant4 Unit Mismatch

**What goes wrong:** Geant4 defaults to mm for position and MeV for energy. The project uses cm and keV. If units are not converted, all LET values will be off by 10x-100x.
**Why it happens:** Geant4 users sometimes change default units in their macro files, so there is no single guaranteed unit system.
**How to avoid:** Accept unit specification as parameters (`position_unit="mm"`, `energy_unit="MeV"`). Apply conversion factors: mm->cm = 0.1, MeV->keV = 1000.
**Warning signs:** LET values outside expected range (0.1-1000 keV/um for therapeutic ions).

### Pitfall 2: Event Boundary Handling

**What goes wrong:** If step-level data is imported but event_id is missing or not properly grouped, the per-event energy sum is wrong.
**Why it happens:** Some Geant4 outputs use sequential step numbering without explicit event_id. Others use a "new event" marker row.
**How to avoid:** Require event_id column. If not present, provide a utility to infer event boundaries (e.g., from step counter resets or position jumps).
**Warning signs:** Suspiciously high or low number of events; event energies that don't match expected spectrum.

### Pitfall 3: LET Outside Table Range

**What goes wrong:** MC events with LET below table minimum or above table maximum get extrapolated incorrectly by `np.interp` (returns endpoint value).
**Why it happens:** The CCE(LET) table covers 0.1-1000 keV/um but some events may fall outside.
**How to avoid:** Log a warning for out-of-range events. `np.interp` clamps to boundary values by default, which is physically reasonable (CCE saturates at high LET, approaches 1 at low LET for fully depleted devices).
**Warning signs:** Cluster of events at exact CCE boundary values.

### Pitfall 4: Empty Events and Zero Energy

**What goes wrong:** Some MC events may deposit zero energy in the SV (particle misses or passes through without interaction). These produce LET=0, which causes log(0) in interpolation.
**Why it happens:** MC transport includes all events, not just those with SV hits.
**How to avoid:** Filter events with edep <= 0 before LET computation. Report the fraction of zero-energy events (useful physics information).
**Warning signs:** Division by zero or -inf in LET calculation.

### Pitfall 5: ROOT File TTree Name Unknown

**What goes wrong:** `uproot.open(path)["wrong_tree_name"]` raises KeyError.
**Why it happens:** Geant4 TTree names are set by the simulation author. There is no universal convention.
**How to avoid:** Provide a `list_trees(filepath)` utility that calls `uproot.open(filepath).keys()` to show available objects. Default tree name should be a reasonable guess but user-overridable.
**Warning signs:** KeyError on ROOT file open.

## Code Examples

### Opening and Inspecting a ROOT File

```python
# Source: uproot official docs (https://uproot.readthedocs.io/en/stable/basic.html)
import uproot

# List all objects in the file
with uproot.open("simulation_output.root") as f:
    print(f.keys())  # e.g., ['Hits;1', 'RunInfo;1']

# Inspect TTree branches
tree = uproot.open("simulation_output.root:Hits")
print(tree.keys())        # branch names
print(tree.typenames())   # branch name -> type mapping
print(tree.num_entries)   # number of entries
```

### Reading TTree Branches as NumPy Arrays

```python
# Source: uproot official docs
import uproot
import numpy as np

tree = uproot.open("simulation_output.root:Hits")

# Read specific branches as numpy
arrays = tree.arrays(["EventID", "PosX", "PosY", "PosZ", "Edep"],
                     library="np")
# arrays["EventID"] -> numpy int32 array
# arrays["Edep"] -> numpy float64 array
```

### Iterating Over Large ROOT Files

```python
# Source: uproot official docs
# For files with millions of entries, iterate in chunks
for batch in uproot.iterate("large_file.root:Hits",
                            ["EventID", "PosX", "Edep"],
                            step_size=10000, library="np"):
    process_batch(batch)
```

### Unit Conversion Pattern

```python
# Geant4 default: mm, MeV -> Project: cm, keV
UNIT_FACTORS = {
    "mm_to_cm": 0.1,
    "cm_to_cm": 1.0,
    "um_to_cm": 1e-4,
    "MeV_to_keV": 1e3,
    "keV_to_keV": 1.0,
    "eV_to_keV": 1e-3,
}

def convert_events(df, pos_unit="mm", energy_unit="MeV"):
    """Convert event DataFrame to project units (cm, keV)."""
    pos_factor = UNIT_FACTORS[f"{pos_unit}_to_cm"]
    energy_factor = UNIT_FACTORS[f"{energy_unit}_to_keV"]

    for col in ["x_cm", "y_cm", "z_cm"]:
        if col in df.columns:
            df[col] = df[col] * pos_factor
    if "edep_keV" in df.columns:
        df["edep_keV"] = df["edep_keV"] * energy_factor

    return df
```

## State of the Art

| Old Approach              | Current Approach      | When Changed       | Impact                                      |
| ------------------------- | --------------------- | ------------------ | ------------------------------------------- |
| PyROOT (C++ dependency)   | uproot (pure Python)  | uproot 4.x (2020+) | No C++ ROOT install needed                  |
| uproot3 (legacy)          | uproot 5.x (current)  | 2022               | API stability, better awkward array support |
| Per-event TCAD simulation | CCE(LET) lookup table | Standard practice  | 10000x speedup for ensemble processing      |
| root_numpy (deprecated)   | uproot                | 2020               | root_numpy no longer maintained             |

**Deprecated/outdated:**

- uproot3: Legacy package, replaced by uproot (5.x). Do not use `import uproot3`.
- root_numpy: Deprecated, removed from scikit-hep. Use uproot instead.

## Open Questions

1. **Geant4 TTree naming from INFN-LNS group**
   - What we know: The group runs Geant4 hadrontherapy simulations. Branch names are user-defined.
   - What's unclear: Exact TTree name, branch names, units, event structure.
   - Recommendation: Make all names configurable with sensible defaults ("Hits" tree, "EventID"/"PosX"/"PosY"/"PosZ"/"Edep" branches). Provide `list_trees()` / `list_branches()` discovery utilities. This blocker is mitigated by configuration flexibility -- it does not block implementation, only final integration testing.

2. **Step-level vs event-total energy deposition**
   - What we know: Geant4 can output per-step (position + edep for each interaction) or per-event (total edep in SV).
   - What's unclear: Which format the INFN-LNS group uses.
   - Recommendation: Support both. Step-level data gets grouped by event_id and summed. Event-total data works directly. Auto-detect based on whether multiple rows share the same event_id.

3. **SV geometry alignment with MC coordinates**
   - What we know: The TCAD mesh has its own coordinate system (x=lateral, y=depth). MC output uses the Geant4 world coordinate system.
   - What's unclear: How the SV is positioned in the Geant4 world.
   - Recommendation: For the CCE(LET) lookup table approach, only the total energy deposited in the SV matters, not the exact position. Position mapping is only needed for the full mesh-resolved approach (MCCP-03 "full path"). Accept an optional coordinate transform parameter for advanced use.

## Sources

### Primary (HIGH confidence)

- uproot official docs (https://uproot.readthedocs.io/en/stable/basic.html) - API patterns, TTree reading, chunked iteration
- Phase 21 `single_particle.py` source code - `load_cce_let_table()`, `ion_track_generation_2d()` API, CCE table JSON format
- Phase 21 `test_single_particle.py` - test patterns for CCE table round-trip

### Secondary (MEDIUM confidence)

- Geant4 hadrontherapy example (https://gitlab.cern.ch/geant4/geant4/-/tree/master/examples/advanced/hadrontherapy/) - typical Geant4 output patterns
- uproot PyPI (https://pypi.org/project/uproot/) - version 4.3.7 (Jan 2026), dependencies
- Geant4 exp_microdosimetry example (https://geant4.web.cern.ch//docs/advanced_examples_doc/example_radioprotection.html) - microdosimetry Geant4 patterns

### Tertiary (LOW confidence)

- INFN-LNS Geant4 TTree naming conventions - unknown, requires sample file from group

## Metadata

**Confidence breakdown:**

- Standard stack: HIGH - uproot and pandas are established, already declared as project dependencies
- Architecture: HIGH - data flow is straightforward (import -> convert -> lookup -> histogram), reuses Phase 21 infrastructure
- Pitfalls: HIGH - well-known issues in MC data coupling (units, event boundaries, LET range)
- ROOT branch naming: LOW - INFN-LNS specific conventions unknown, mitigated by configurable design

**Research date:** 2026-03-31
**Valid until:** 2026-04-30 (stable domain, no fast-moving dependencies)

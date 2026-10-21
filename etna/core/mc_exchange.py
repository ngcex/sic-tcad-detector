"""File exchange with Allpix Squared (and Geant4-style hit lists). No Allpix/ROOT needed.

1. Hit lists (energy deposits) both ways, in the Allpix Squared ``DepositionReader``
   CSV format (module README and ``DepositionReaderModule::read_csv``, Allpix Squared
   master, read 2026-10-21):

       Event: <N>
       <PDG>, <time>, <energy>, <x>, <y>, <z>, <detector>, <track_id>, <parent_id>

   ``#`` lines and blank lines are ignored; an event header line starts with ``E``;
   events are numbered from 0; units are set in the Allpix config (``unit_length``,
   ``unit_time``, ``unit_energy``; framework defaults mm / ns / MeV); parent id 0 marks
   the primary and a parent must appear before it is used.

   - :func:`write_deposition_csv` writes ETNA-sampled tracks (``sample_track_events(...,
     max_hits=...)``) so Allpix can propagate them in its own fields;
   - :func:`read_deposition_csv` reads such a file (from Allpix users, or a Geant4
     standalone run written in this format) into the step-level DataFrame of
     ``etna.core.mc_coupling`` (event_id, x_cm, y_cm, z_cm, edep_keV, ...), i.e.
     straight into ``process_mc_ensemble`` and the CCE chain.

2. Field maps in the Allpix Squared legacy INIT format, exactly as
   ``FieldParser::parse_init_file`` / ``FieldWriter::write_init_file``
   (``src/tools/field_parser.h``) read and write them:

       <header line>
       <units> <norm>
       ##TURN## ##TILT## 1.0
       0.0 0.0 0.0
       <thickness_um> <size_x_um> <size_y_um> 0.0 0.0 0.0 0.0 <nx> <ny> <nz> 0.0
       i j k v1 [v2 v3]       (1-based; x outer loop, y, z inner)

   The parser treats the value grid as CELLS: z index k covers
   [z_lo + (k-1) d/nz, z_lo + k d/nz] of the field depth d, and z_lo is the BACK of
   the sensor (Allpix local +z points to the implant side). ElectricFieldReader reads
   V/cm (``file_units`` default), DopingProfileReader 1/cm^3 (positive = n-type),
   WeightingPotentialReader a dimensionless scalar.

   :func:`export_allpix_fields` converts an ETNA ``run_field`` result: the potential is
   interpolated onto the cell centres and differentiated (E = -grad psi), so signs are
   physical; ETNA depth 0 (the p+ anode / implant side) maps to Allpix z = +t/2.

Not verified: no Allpix Squared installation was run against these files (the
writer/reader are checked against each other and against the quoted parser layout,
including an asymmetric-grid index test). The global-coordinate placement of the SV
in an Allpix geometry is the user's (see ``write_deposition_csv``).
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

LENGTH_TO_CM = {"nm": 1e-7, "um": 1e-4, "mm": 0.1, "cm": 1.0, "m": 100.0}
ENERGY_TO_KEV = {"eV": 1e-3, "keV": 1.0, "MeV": 1e3, "GeV": 1e6}
TIME_TO_NS = {"ps": 1e-3, "ns": 1.0, "us": 1e3, "s": 1e9}
PDG_ELECTRON = 11


# ------------------------------------------------------------------ hit lists ----


def tracks_to_hits(events, beam_from_implant_side=True):
    """Hit table (dict of arrays) from ``sample_track_events(..., max_hits>0)``.

    Coordinates (um) in an SV-centred frame: x, y as sampled, z' = t/2 - z when
    ``beam_from_implant_side`` (the beam enters at the Allpix +z / implant face and
    travels to -z), else z' = z - t/2. Track id 1 = primary (parent 0), delta rays
    2, 3, ... (parent 1; parent 0 if the primary left no deposit in the SV).
    """
    if events.hits is None:
        raise ValueError("events carry no hits: sample with max_hits > 0")
    t = events.meta["sv_um"][2]
    pdg_primary = {"proton": 2212, "alpha": 1000020040, "c12": 1000060120}[
        events.meta["particle"]
    ]
    rows = {
        k: [] for k in ("event", "pdg", "x", "y", "z", "edep_keV", "track", "parent")
    }
    nh = events.n_hits[:, 0]
    for i in np.nonzero(nh)[0]:
        h = events.hits[i, : nh[i]]
        trk = h[:, 4].astype(int)
        has_primary = np.any(trk == 1)
        z = (t / 2.0 - h[:, 2]) if beam_from_implant_side else (h[:, 2] - t / 2.0)
        order = np.argsort(trk, kind="stable")  # primary first
        for j in order:
            rows["event"].append(i)
            rows["pdg"].append(pdg_primary if trk[j] == 1 else PDG_ELECTRON)
            rows["x"].append(h[j, 0])
            rows["y"].append(h[j, 1])
            rows["z"].append(z[j])
            rows["edep_keV"].append(h[j, 3])
            rows["track"].append(trk[j])
            rows["parent"].append(0 if trk[j] == 1 else (1 if has_primary else 0))
    return {k: np.asarray(v) for k, v in rows.items()}


def write_deposition_csv(
    path,
    events,
    detector="sic_sv",
    unit_length="um",
    unit_energy="keV",
    beam_from_implant_side=True,
    origin_um=(0.0, 0.0, 0.0),
):
    """Write sampled tracks as an Allpix Squared DepositionReader CSV file.

    Every event of ``events`` gets an ``Event: <index>`` header (also events without
    deposits, so the Allpix event numbering stays aligned with ETNA's). Positions are
    the SV-centred frame of :func:`tracks_to_hits` shifted by ``origin_um`` (place the
    Allpix detector so its sensor centre is at ``origin_um``, unrotated; Allpix local
    +z towards the implants). Time column = 0 (the sampler has no timing).

    Returns the matching ``[DepositionReader]`` config snippet (str).
    """
    if unit_length not in LENGTH_TO_CM or unit_energy not in ENERGY_TO_KEV:
        raise ValueError("unsupported unit_length / unit_energy")
    hits = tracks_to_hits(events, beam_from_implant_side)
    lf = 1e-4 / LENGTH_TO_CM[unit_length]  # um -> unit
    ef = 1.0 / ENERGY_TO_KEV[unit_energy]
    ox, oy, oz = (float(v) for v in origin_um)
    lines = [
        "# ETNA mc_transport hits (etna.core.mc_exchange.write_deposition_csv)",
        f"# particle={events.meta['particle']} sv_um={events.meta['sv_um']} "
        f"tcut_keV={events.meta['tcut_keV']:.4g} seed={events.meta['seed']}",
        f"# units: length {unit_length}, time ns, energy {unit_energy}",
    ]
    ev = hits["event"]
    starts = np.searchsorted(ev, np.arange(events.n_events + 1))
    for i in range(events.n_events):
        lines.append(f"Event: {i}")
        for j in range(starts[i], starts[i + 1]):
            lines.append(
                f"{hits['pdg'][j]}, 0, {hits['edep_keV'][j] * ef:.9e}, "
                f"{(hits['x'][j] + ox) * lf:.9e}, {(hits['y'][j] + oy) * lf:.9e}, "
                f"{(hits['z'][j] + oz) * lf:.9e}, {detector}, {hits['track'][j]}, "
                f"{hits['parent'][j]}"
            )
    Path(path).write_text("\n".join(lines) + "\n", encoding="utf-8")  # EOF on new line
    return (
        "[DepositionReader]\n"
        'model = "csv"\n'
        f'file_name = "{Path(path).name}"\n'
        f'unit_length = "{unit_length}"\n'
        'unit_time = "ns"\n'
        f'unit_energy = "{unit_energy}"\n'
    )


def read_deposition_csv(
    path,
    unit_length="mm",
    unit_energy="MeV",
    unit_time="ns",
    detector=None,
    assign_timestamps=True,
    create_mcparticles=True,
):
    """Read an Allpix Squared DepositionReader CSV into an ETNA step DataFrame.

    Units default to the Allpix framework defaults (mm, MeV, ns) like the module;
    ``assign_timestamps`` / ``create_mcparticles`` False mean the time / track / parent
    columns are absent (as in Allpix). ``detector`` keeps only that volume name.
    Columns: event_id, pdg, time_ns, x_cm, y_cm, z_cm, edep_keV, detector, track_id,
    parent_id (``track_id`` == 1 is NOT assumed to be the primary: parent_id == 0 is).
    """
    import pandas as pd

    for u, table in (
        (unit_length, LENGTH_TO_CM),
        (unit_energy, ENERGY_TO_KEV),
        (unit_time, TIME_TO_NS),
    ):
        if u not in table:
            raise ValueError(f"unsupported unit {u!r}")
    rows = []
    event = None
    for raw in Path(path).read_text(encoding="utf-8").splitlines():
        line = raw.strip()
        if not line or line.startswith("#"):
            continue
        if line[0] == "E":
            parts = line.split()
            event = int(parts[1])
            continue
        if event is None:
            raise ValueError("deposit line before the first 'Event:' header")
        f = [c.strip() for c in line.split(",")]
        k = 0
        pdg = int(f[k])
        k += 1
        t = 0.0
        if assign_timestamps:
            t = float(f[k])
            k += 1
        e = float(f[k])
        x, y, z = float(f[k + 1]), float(f[k + 2]), float(f[k + 3])
        vol = f[k + 4]
        trk = par = 0
        if create_mcparticles:
            trk, par = int(f[k + 5]), int(f[k + 6])
        rows.append((event, pdg, t, x, y, z, e, vol, trk, par))
    df = pd.DataFrame(
        rows,
        columns=[
            "event_id",
            "pdg",
            "time_ns",
            "x_cm",
            "y_cm",
            "z_cm",
            "edep_keV",
            "detector",
            "track_id",
            "parent_id",
        ],
    )
    if detector is not None:
        df = df[df["detector"] == detector].reset_index(drop=True)
    lf = LENGTH_TO_CM[unit_length]
    for c in ("x_cm", "y_cm", "z_cm"):
        df[c] = df[c] * lf
    df["edep_keV"] = df["edep_keV"] * ENERGY_TO_KEV[unit_energy]
    df["time_ns"] = df["time_ns"] * TIME_TO_NS[unit_time]
    return df


def tracks_to_dataframe(events):
    """Event-level DataFrame of sampled tracks for ``process_mc_ensemble``: one row per
    event with edep_keV, the SV-entry position (cm), ``ekin_MeV`` = kinetic energy at SV
    entry (for the kappa(E) lookup) and ``track_id`` = 1."""
    import pandas as pd

    ent = np.nan_to_num(events.entry_um, nan=0.0) * 1e-4
    return pd.DataFrame(
        {
            "event_id": np.arange(events.n_events),
            "x_cm": ent[:, 0],
            "y_cm": ent[:, 1],
            "z_cm": ent[:, 2],
            "edep_keV": events.edep_keV,
            "ekin_MeV": events.E_sv_entry_MeV,
            "track_id": np.ones(events.n_events, dtype=int),
        }
    )


# ------------------------------------------------------------------ INIT fields ----


def write_init(path, values, size_um, units, header="ETNA field", norm=1.0):
    """Write a regular field grid in the Allpix INIT format.

    values: (nx, ny, nz) scalar or (nx, ny, nz, 3) vector grid, index [i, j, k] = cell
    (x_i, y_j, z_k), z from the back (k = 0) to the implant side. size_um =
    (size_x, size_y, thickness). ``units``: e.g. "V/cm", "/cm/cm/cm"; "" writes the
    Allpix placeholder "internal" (dimensionless, weighting potential).
    """
    v = np.asarray(values, dtype=float)
    if v.ndim == 3:
        v = v[..., None]
    if v.ndim != 4 or v.shape[3] not in (1, 3):
        raise ValueError("values must be (nx, ny, nz) or (nx, ny, nz, 3)")
    if not np.all(np.isfinite(v)):
        raise ValueError("field values must be finite")
    nx, ny, nz, nc = v.shape
    sx, sy, th = (float(s) for s in size_um)
    if "\n" in header:
        raise ValueError("header must be one line")
    out = [
        header,
        f"{units if units else 'internal'} {norm:g}",
        "##TURN## ##TILT## 1.0",
        "0.0 0.0 0.0",
        f"{th:.10g} {sx:.10g} {sy:.10g} 0.0 0.0 0.0 0.0 {nx} {ny} {nz} 0.0",
    ]
    idx = np.indices((nx, ny, nz)).reshape(3, -1).T + 1
    flat = v.reshape(-1, nc)
    for (i, j, k), row in zip(idx, flat):
        out.append(f"{i} {j} {k} " + " ".join(f"{x:.10e}" for x in row))
    Path(path).write_text("\n".join(out) + "\n", encoding="utf-8")


def read_init(path):
    """Read an INIT file the way ``FieldParser::parse_init_file`` does.

    Returns dict: header, units, norm, size_um (x, y, thickness), shape (nx, ny, nz),
    values (nx, ny, nz, n_components).
    """
    with open(path, encoding="utf-8") as fh:
        header = fh.readline().rstrip("\n")
        tok = fh.read().split()
    units = tok[0]
    try:
        norm = float(tok[1])
    except ValueError:
        norm = 1.0
    # tok[2:5] turn/tilt, tok[5:8] magnetic field
    th, sx, sy = float(tok[8]), float(tok[9]), float(tok[10])
    nx, ny, nz = int(tok[15]), int(tok[16]), int(tok[17])
    data = np.asarray(tok[19:], dtype=float)
    nvert = nx * ny * nz
    nc = data.size // nvert - 3
    if nc not in (1, 3) or data.size != nvert * (3 + nc):
        raise ValueError("unexpected INIT data block size")
    data = data.reshape(nvert, 3 + nc)
    vals = np.full((nx, ny, nz, nc), np.nan)
    ii = data[:, :3].astype(int) - 1
    if ii.min() < 0 or np.any(ii.max(axis=0) >= (nx, ny, nz)):
        raise ValueError("INIT index out of range")
    vals[ii[:, 0], ii[:, 1], ii[:, 2]] = data[:, 3:]
    return {
        "header": header,
        "units": units,
        "norm": norm,
        "size_um": (sx, sy, th),
        "shape": (nx, ny, nz),
        "values": vals,
    }


def _cell_depths_um(thickness_um, nz):
    """ETNA depth (um, 0 at the implant side) of the Allpix z-cell centres k = 0..nz-1
    (k = 0 is the back)."""
    zc = (np.arange(nz) + 0.5) / nz * thickness_um  # from the back
    return thickness_um - zc


def export_allpix_fields(
    result,
    outdir,
    nz=200,
    nx=None,
    pitch_um=(100.0, 100.0),
    depth_range_um=None,
    prefix="etna",
):
    """Write ``<prefix>_efield.init`` (V/cm), ``<prefix>_doping.init`` (1/cm^3) and
    ``<prefix>_weighting.init`` (planar pad) from a ``run_field`` SimResult.

    1D results (mesh x = depth): a laterally uniform map (nx = ny = 1, Allpix uses the
    single bin everywhere), pitch_um sets the declared lateral size. 2D results (mesh
    x = lateral, y = depth): the potential is interpolated onto an nx x nz grid over
    the mesh's lateral range (mirrored to a full plane if it starts at 0, i.e. a
    half-cell model), ny = 1. The weighting potential is the planar-pad 1D solution
    phi_w = z/t (1 at the implant side, 0 at the back): exact for a pad much wider
    than thick, an approximation near the pad edge.

    Returns dict: paths, thickness_um, and an Allpix config snippet.
    """
    mesh = getattr(result, "mesh", None)
    if mesh is None or "Potential" not in mesh.node_values:
        raise ValueError("result has no mesh Potential (use run_field)")
    out = Path(outdir)
    out.mkdir(parents=True, exist_ok=True)
    psi = np.asarray(mesh.node_values["Potential"], dtype=float)
    dop = np.asarray(mesh.node_values.get("NetDoping", np.zeros_like(psi)), dtype=float)
    two_d = (
        mesh.y_coords is not None
        and np.size(mesh.y_coords) == psi.size
        and (np.ptp(np.asarray(mesh.y_coords)) > 0)
    )
    if not two_d:
        depth = np.asarray(mesh.x_coords, dtype=float) * 1e4
        o = np.argsort(depth)
        depth, psi_s, dop_s = depth[o], psi[o], dop[o]
        d0, d1 = (depth[0], depth[-1]) if depth_range_um is None else depth_range_um
        th = float(d1 - d0)
        zdep = d0 + _cell_depths_um(th, nz)
        # E_z(Allpix) = -dpsi/dz = +dpsi/d(depth): central differences on the cell grid
        h = th / nz
        Ez = (
            np.interp(zdep + 0.5 * h, depth, psi_s)
            - np.interp(zdep - 0.5 * h, depth, psi_s)
        ) / (h * 1e-4)
        ef = np.zeros((1, 1, nz, 3))
        ef[0, 0, :, 2] = Ez
        dp = np.interp(zdep, depth, dop_s).reshape(1, 1, nz)
        size = (float(pitch_um[0]), float(pitch_um[1]), th)
        mapping = "SENSOR"
    else:
        from scipy.interpolate import LinearNDInterpolator

        xl = np.asarray(mesh.x_coords, dtype=float) * 1e4
        yd = np.asarray(mesh.y_coords, dtype=float) * 1e4
        d0, d1 = (yd.min(), yd.max()) if depth_range_um is None else depth_range_um
        th = float(d1 - d0)
        mirror = abs(xl.min()) < 1e-6
        x_hi = xl.max()
        x_lo = -x_hi if mirror else xl.min()
        nx = 101 if nx is None else int(nx)
        xc = x_lo + (np.arange(nx) + 0.5) / nx * (x_hi - x_lo)
        zdep = d0 + _cell_depths_um(th, nz)
        ip = LinearNDInterpolator(np.column_stack([xl, yd]), psi)
        idop = LinearNDInterpolator(np.column_stack([xl, yd]), dop)
        X, Z = np.meshgrid(np.abs(xc) if mirror else xc, zdep, indexing="ij")
        hx = (x_hi - x_lo) / nx
        hz = th / nz

        def _p(f, a, b):
            a = np.clip(a, xl.min(), xl.max())
            b = np.clip(b, yd.min(), yd.max())
            v = f(a, b)
            if np.any(~np.isfinite(v)):
                from scipy.interpolate import NearestNDInterpolator

                nn = NearestNDInterpolator(np.column_stack([xl, yd]), f.values.ravel())
                v = np.where(np.isfinite(v), v, nn(a, b))
            return v

        Ex = -(_p(ip, X + 0.5 * hx, Z) - _p(ip, X - 0.5 * hx, Z)) / (hx * 1e-4)
        if mirror:
            Ex = Ex * np.sign(xc)[:, None]
        Ez = (_p(ip, X, Z + 0.5 * hz) - _p(ip, X, Z - 0.5 * hz)) / (hz * 1e-4)
        ef = np.zeros((nx, 1, nz, 3))
        ef[:, 0, :, 0] = Ex
        ef[:, 0, :, 2] = Ez
        dp = _p(idop, X, Z).reshape(nx, 1, nz)
        size = (float(x_hi - x_lo), float(pitch_um[1]), th)
        mapping = "PIXEL_FULL"
    wp = ((np.arange(nz) + 0.5) / nz).reshape(1, 1, nz)
    paths = {
        "efield": out / f"{prefix}_efield.init",
        "doping": out / f"{prefix}_doping.init",
        "weighting": out / f"{prefix}_weighting.init",
    }
    cfg = result.config
    hdr = f"ETNA run_field export (bias {result.metadata.get('bias_V', 'n/a')} V)"
    write_init(paths["efield"], ef, size, "V/cm", header=hdr)
    write_init(paths["doping"], dp, size, "/cm/cm/cm", header=hdr)
    write_init(
        paths["weighting"], wp, (size[0], size[1], th), "", header=hdr + " planar pad"
    )
    snippet = (
        "[ElectricFieldReader]\n"
        'model = "mesh"\n'
        f'file_name = "{paths["efield"].name}"\n'
        f'field_mapping = "{mapping}"\n\n'
        "[DopingProfileReader]\n"
        'model = "mesh"\n'
        f'file_name = "{paths["doping"].name}"\n'
        f'field_mapping = "{mapping}"\n\n'
        "[WeightingPotentialReader]\n"
        'model = "mesh"\n'
        f'file_name = "{paths["weighting"].name}"\n'
        'field_mapping = "PIXEL_FULL"\n'
    )
    return {
        "paths": {k: str(v) for k, v in paths.items()},
        "thickness_um": th,
        "size_um": size,
        "config_snippet": snippet,
        "config_name": type(cfg).__name__,
    }

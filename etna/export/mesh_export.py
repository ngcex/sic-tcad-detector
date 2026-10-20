"""Per-node profile export of an ``etna.run_field`` result: CSV, legacy-ASCII VTK, VTU.

Quantities per node: NetDoping (cm^-3, = N_D - N_A with IONIZED acceptors, the devsim
NetDoping model), Potential (V), ElectricField magnitude (V/cm) and, when the result
carries them (``run_field`` does since 2026-10-20), the carrier densities Electrons and
Holes (cm^-3). devsim path: the node solutions of the solved device. 1D batched-engine
path: a dark Gummel solve on the Poisson state, within 2e-9 relative of devsim at every
node (tests/test_field_engine_1d.py). 2D engine path:
the engine dark state (``etna.core.engine2d.dark_state_2d``); measured against devsim on
the default stack, hw 25 um, -20 V: same values to |dlog10| < 1e-13 wherever the
density exceeds 1e10 cm^-3; in the depleted region both are < 1e-8 cm^-3 (devsim
~1e-13, the engine underflows towards 0), i.e. physically zero in both. Older results
without carriers export the three base columns only. Coordinates are written in um. 1D: x = depth from the anode. 2D: x = lateral coordinate, y = depth from the
anode (etna.core.device2d convention). NaN values are written as ``nan``.

VTK (hand-written, no dependency): 1D -> RECTILINEAR_GRID (N x 1 x 1). 2D -> a
RECTILINEAR_GRID when the nodes form a tensor-product grid (the 2D mesh is; the node
order from the solver is not grid order, so every array is permuted with the same
x-fastest lexsort), otherwise POLYDATA vertices.

VTU (XML UnstructuredGrid, ASCII, hand-written): 1D -> VTK_LINE cells along x; 2D
tensor grid -> VTK_QUAD cells; otherwise VTK_VERTEX cells. Same arrays, units and
label (FieldData) as the legacy file. Readable by ParaView / VisIt / meshio.
"""

from __future__ import annotations

import csv
import io
from xml.sax.saxutils import escape, quoteattr

import numpy as np

from etna.export.labels import VALIDATION_LABEL

_FIELDS = (
    ("NetDoping", "NetDoping_cm-3"),
    ("Potential", "Potential_V"),
    ("ElectricField", "ElectricField_V_per_cm"),
)
# exported when present in MeshData.node_values (appended after the base columns)
_OPTIONAL_FIELDS = (
    ("Electrons", "Electrons_cm-3"),
    ("Holes", "Holes_cm-3"),
)


def _arrays(result):
    if result.mesh is None:
        raise ValueError("result has no mesh (use the SimResult of etna.run_field)")
    m = result.mesh
    x = np.asarray(m.x_coords, dtype=float) * 1e4  # cm -> um
    y = None if m.y_coords is None else np.asarray(m.y_coords, dtype=float) * 1e4
    vals = {k: np.asarray(m.node_values[k], dtype=float) for k, _ in _FIELDS}
    for k, _ in _OPTIONAL_FIELDS:
        if k in m.node_values:
            vals[k] = np.asarray(m.node_values[k], dtype=float)
    for k, v in vals.items():
        if v.shape != x.shape:
            raise ValueError(f"{k} has {v.size} values for {x.size} nodes")
    return x, y, vals


def _fields(vals) -> "list[tuple[str, str]]":
    return list(_FIELDS) + [(k, h) for k, h in _OPTIONAL_FIELDS if k in vals]


def _meta_lines(result) -> "list[str]":
    md = result.metadata
    two_d = result.mesh is not None and result.mesh.y_coords is not None
    return [
        VALIDATION_LABEL,
        f"solver: {md.get('solver', '')}; validated_domain: {md.get('validated_domain')}; "
        f"domain_note: {md.get('domain_note') or ''}"
        + (f"; fallback: {md['fallback_reason']}" if md.get("fallback_reason") else ""),
        f"bias_V: {md.get('bias_V')} (negative = reverse bias)",
        (
            "coordinates um: x = lateral, y = depth from the anode"
            if two_d
            else "coordinates um: x = depth from the anode"
        ),
        "NetDoping = N_D - N_A (ionized acceptors), cm^-3; Potential V; "
        "ElectricField = |E| V/cm; Electrons / Holes cm^-3 (when present); "
        "NaN written as nan",
    ]


def profile_to_csv(result) -> str:
    """CSV with '#' comment header (label, solver, validity, units) then one row per node."""
    x, y, vals = _arrays(result)
    buf = io.StringIO()
    for line in _meta_lines(result):
        buf.write(f"# {line}\n")
    w = csv.writer(buf)
    fields = _fields(vals)
    head = ["x_um"] + (["y_um"] if y is not None else []) + [h for _, h in fields]
    w.writerow(head)
    cols = [x] + ([y] if y is not None else []) + [vals[k] for k, _ in fields]
    for row in zip(*cols):
        w.writerow([repr(float(v)) if np.isfinite(v) else "nan" for v in row])
    return buf.getvalue()


def _num(v) -> str:
    return repr(float(v)) if np.isfinite(v) else "nan"


def _tensor_order(x, y):
    """Permutation to x-fastest grid order and the unique axes, or None if not a grid."""
    ux, uy = np.unique(x), np.unique(y)
    if ux.size * uy.size != x.size:
        return None
    if np.unique(np.stack([x, y], axis=1), axis=0).shape[0] != x.size:
        return None
    order = np.lexsort((x, y))  # primary key y, secondary x -> x varies fastest
    if not (
        np.array_equal(x[order], np.tile(ux, uy.size))
        and np.array_equal(y[order], np.repeat(uy, ux.size))
    ):
        return None
    return order, ux, uy


def profile_to_vtk(result) -> str:
    """Legacy ASCII VTK (version 3.0) string of the node profile (1D or 2D)."""
    x, y, vals = _arrays(result)
    md = result.metadata
    title = (
        f"{VALIDATION_LABEL} | etna run_field bias {md.get('bias_V')} V, "
        f"solver {md.get('solver', '')}, coords um"
    ).replace("\n", " ")[:255]
    out = ["# vtk DataFile Version 3.0", title, "ASCII"]
    if y is None:
        o = np.argsort(x, kind="stable")
        out += ["DATASET RECTILINEAR_GRID", f"DIMENSIONS {x.size} 1 1"]
        out += [f"X_COORDINATES {x.size} double", " ".join(_num(v) for v in x[o])]
        out += ["Y_COORDINATES 1 double", "0", "Z_COORDINATES 1 double", "0"]
        n = x.size
    else:
        t = _tensor_order(x, y)
        if t is not None:
            o, ux, uy = t
            out += ["DATASET RECTILINEAR_GRID", f"DIMENSIONS {ux.size} {uy.size} 1"]
            out += [f"X_COORDINATES {ux.size} double", " ".join(_num(v) for v in ux)]
            out += [f"Y_COORDINATES {uy.size} double", " ".join(_num(v) for v in uy)]
            out += ["Z_COORDINATES 1 double", "0"]
        else:
            o = np.arange(x.size)
            out += ["DATASET POLYDATA", f"POINTS {x.size} double"]
            out += [f"{_num(a)} {_num(b)} 0" for a, b in zip(x, y)]
            out += [f"VERTICES {x.size} {2 * x.size}"] + [
                f"1 {i}" for i in range(x.size)
            ]
        n = x.size
    validated = md.get("validated_domain")
    out += [
        "FIELD FieldData 2",
        "bias_V 1 1 double",
        _num(float("nan") if md.get("bias_V") is None else md["bias_V"]),
        "validated_domain 1 1 int",
        "1" if validated else "0",
        f"POINT_DATA {n}",
    ]
    for key, name in _fields(vals):
        v = vals[key][o]
        out += [f"SCALARS {name} double 1", "LOOKUP_TABLE default"]
        out += [" ".join(_num(a) for a in v[i : i + 8]) for i in range(0, v.size, 8)]
    return "\n".join(out) + "\n"


# VTK cell type ids (vtkCellType.h)
VTK_VERTEX, VTK_LINE, VTK_QUAD = 1, 3, 9


def _cells(x, y):
    """(order, connectivity rows, cell type) on the reordered points."""
    if y is None:
        o = np.argsort(x, kind="stable")
        conn = [(i, i + 1) for i in range(x.size - 1)]
        return o, conn, VTK_LINE
    t = _tensor_order(x, y)
    if t is None:
        return np.arange(x.size), [(i,) for i in range(x.size)], VTK_VERTEX
    o, ux, uy = t
    nx = ux.size
    conn = [
        (j * nx + i, j * nx + i + 1, (j + 1) * nx + i + 1, (j + 1) * nx + i)
        for j in range(uy.size - 1)
        for i in range(nx - 1)
    ]
    return o, conn, VTK_QUAD


def _rows(values, per_line=6) -> str:
    vals = list(values)
    return "\n".join(
        " ".join(vals[i : i + per_line]) for i in range(0, len(vals), per_line)
    )


def profile_to_vtu(result) -> str:
    """VTK XML UnstructuredGrid (ASCII) string of the node profile (1D or 2D).

    Points in um (z = 0); point data as in ``profile_to_vtk``; FieldData carries
    bias_V and validated_domain; the validation label, solver and units are in the XML
    comment after the root element.
    """
    x, y, vals = _arrays(result)
    md = result.metadata
    o, conn, ctype = _cells(x, y)
    yy = np.zeros_like(x) if y is None else y
    pts = _rows((f"{_num(a)} {_num(b)} 0" for a, b in zip(x[o], yy[o])), 3)
    offsets, acc = [], 0
    for c in conn:
        acc += len(c)
        offsets.append(acc)
    validated = md.get("validated_domain")
    bias = float("nan") if md.get("bias_V") is None else md["bias_V"]
    meta = " | ".join(_meta_lines(result)).replace("--", "- -")
    out = [
        '<?xml version="1.0"?>',
        '<VTKFile type="UnstructuredGrid" version="1.0" byte_order="LittleEndian" '
        'header_type="UInt64">',
        "<!-- " + escape(meta) + " -->",
        "<UnstructuredGrid>",
        "<FieldData>",
        '<DataArray type="Float64" Name="bias_V" NumberOfTuples="1" format="ascii">'
        f"{_num(bias)}</DataArray>",
        '<DataArray type="Int32" Name="validated_domain" NumberOfTuples="1" '
        f'format="ascii">{1 if validated else 0}</DataArray>',
        "</FieldData>",
        f'<Piece NumberOfPoints="{x.size}" NumberOfCells="{len(conn)}">',
        "<PointData>",
    ]
    for key, name in _fields(vals):
        out += [
            f'<DataArray type="Float64" Name={quoteattr(name)} format="ascii">',
            _rows(_num(a) for a in vals[key][o]),
            "</DataArray>",
        ]
    out += [
        "</PointData>",
        "<Points>",
        '<DataArray type="Float64" NumberOfComponents="3" format="ascii">',
        pts,
        "</DataArray>",
        "</Points>",
        "<Cells>",
        '<DataArray type="Int64" Name="connectivity" format="ascii">',
        _rows((" ".join(str(i) for i in c) for c in conn), 4),
        "</DataArray>",
        '<DataArray type="Int64" Name="offsets" format="ascii">',
        _rows((str(v) for v in offsets), 12),
        "</DataArray>",
        '<DataArray type="UInt8" Name="types" format="ascii">',
        _rows((str(ctype) for _ in conn), 24),
        "</DataArray>",
        "</Cells>",
        "</Piece>",
        "</UnstructuredGrid>",
        "</VTKFile>",
    ]
    return "\n".join(out) + "\n"

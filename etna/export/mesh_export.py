"""Per-node profile export of an ``etna.run_field`` result: CSV and legacy-ASCII VTK.

Quantities per node: NetDoping (cm^-3, = N_D - N_A with IONIZED acceptors, the devsim
NetDoping model), Potential (V), ElectricField magnitude (V/cm). Coordinates are written
in um. 1D: x = depth from the anode. 2D: x = lateral coordinate, y = depth from the
anode (etna.core.device2d convention). NaN values are written as ``nan``.

VTK (hand-written, no dependency): 1D -> RECTILINEAR_GRID (N x 1 x 1). 2D -> a
RECTILINEAR_GRID when the nodes form a tensor-product grid (the 2D mesh is; the node
order from the solver is not grid order, so every array is permuted with the same
x-fastest lexsort), otherwise POLYDATA vertices.
"""

from __future__ import annotations

import csv
import io

import numpy as np

from etna.export.labels import VALIDATION_LABEL

_FIELDS = (
    ("NetDoping", "NetDoping_cm-3"),
    ("Potential", "Potential_V"),
    ("ElectricField", "ElectricField_V_per_cm"),
)


def _arrays(result):
    if result.mesh is None:
        raise ValueError("result has no mesh (use the SimResult of etna.run_field)")
    m = result.mesh
    x = np.asarray(m.x_coords, dtype=float) * 1e4  # cm -> um
    y = None if m.y_coords is None else np.asarray(m.y_coords, dtype=float) * 1e4
    vals = {k: np.asarray(m.node_values[k], dtype=float) for k, _ in _FIELDS}
    for k, v in vals.items():
        if v.shape != x.shape:
            raise ValueError(f"{k} has {v.size} values for {x.size} nodes")
    return x, y, vals


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
        "ElectricField = |E| V/cm; NaN written as nan",
    ]


def profile_to_csv(result) -> str:
    """CSV with '#' comment header (label, solver, validity, units) then one row per node."""
    x, y, vals = _arrays(result)
    buf = io.StringIO()
    for line in _meta_lines(result):
        buf.write(f"# {line}\n")
    w = csv.writer(buf)
    head = ["x_um"] + (["y_um"] if y is not None else []) + [h for _, h in _FIELDS]
    w.writerow(head)
    cols = [x] + ([y] if y is not None else []) + [vals[k] for k, _ in _FIELDS]
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
    for key, name in _FIELDS:
        v = vals[key][o]
        out += [f"SCALARS {name} double 1", "LOOKUP_TABLE default"]
        out += [" ".join(_num(a) for a in v[i : i + 8]) for i in range(0, v.size, 8)]
    return "\n".join(out) + "\n"

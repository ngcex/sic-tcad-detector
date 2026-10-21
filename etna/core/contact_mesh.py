"""Depth-mesh grading of the contact-on-epi builders (cathode contact cell and junction).

Why (issue #4 follow-up, lab/mesh_contact/notes/LOG.md). An ohmic contact on the n- epi
is a perfect sink for holes (the minority carrier there). Holes generated in its boundary
layer are lost: the layer is the undepleted epi next to the contact at partial depletion,
or the drift-diffusion layer of thickness ~ kT / (q E) (13-26 nm at 1-2e4 V/cm) once the
depletion reaches the contact. The box scheme lumps the contact cell's generation onto the
contact node, where the contact equation replaces the continuity equation, so the CCE
converges at FIRST order in the contact-cell size and a coarse cathode cell OVERESTIMATES
the CCE: +6.2e-3 at 10 V and +7.9e-4 at 50 V on the 100 nm cell of the
alternative-structure builders (mesa centre, against the converged 0.9630651 / 0.9990906).
The next error at partial depletion is the 10 nm junction spacing of those builders: the
CCE at 10 V sits ~0.5 um from reach-through and moves with the depletion edge.

``"graded"`` (default since 2026-10-21): the cathode cell is 0.1 nm, the cells grow by
1.25x away from the contact up to 50 nm and stay 50 nm over the last 1 um of the epi; the
junction line has 1 nm spacing (as device2d / the 1D mesh). Mesa centre CCE error against
the converged value: +5.4e-5 at 10 V, -2.9e-6 at 50 V (legacy: +6.2e-3 / +7.9e-4).
``"legacy"``: the meshes before 2026-10-21 (100 nm cathode cell, 10 nm junction), for
reproducing earlier numbers.
A dict overrides single keys of ``DEPTH_MESH_GRADED`` (convergence studies).

The default n+ stack (DEFAULT_STACK / DeviceConfig()) has no such contact: its cathode
sits on the 1e19 n+ layer, the generation is confined to the epi and the n-/n+ high-low
junction reflects holes. The 1D engine there moves by < 2e-10 when its contact cells are
bisected (lab/mesh_contact/bench/default_1d.py), so the default meshes are not changed.
"""

from __future__ import annotations

DEPTH_MESH_LEGACY = "legacy"
DEPTH_MESH_DEFAULT = "graded"

# The legacy depth mesh of etna.core.alternative_structures (cm).
LEGACY_CATHODE_PS = 1e-5
LEGACY_JUNCTION_PS = 1e-6

# Default grading (convergence study: lab/mesh_contact/notes/LOG.md section 3).
# first_cell_cm: the contact cell; ratio: geometric growth away from the contact until
# the cell reaches near_ps_cm; then uniform near_ps_cm cells up to near_depth_cm from the
# contact (the undepleted epi next to the cathode at partial depletion); beyond that the
# builder's own lines take over (devsim grades from near_ps_cm to the epi spacing).
# junction_ps_cm: spacing of the p+/n- junction line.
DEPTH_MESH_GRADED = {
    "first_cell_cm": 1e-8,
    "ratio": 1.25,
    "near_ps_cm": 5e-6,
    "near_depth_cm": 1e-4,
    "junction_ps_cm": 1e-7,
}


def resolve_depth_mesh(spec=DEPTH_MESH_DEFAULT):
    """None for the legacy mesh, else the grading dict (defaults overridden by ``spec``).

    ``spec``: "graded" (default), "legacy", or a dict with any of the keys of
    ``DEPTH_MESH_GRADED``.
    """
    if spec is None or spec == DEPTH_MESH_LEGACY:
        return None
    if spec == DEPTH_MESH_DEFAULT:
        return dict(DEPTH_MESH_GRADED)
    if isinstance(spec, dict):
        unknown = set(spec) - set(DEPTH_MESH_GRADED)
        if unknown:
            raise ValueError(f"depth_mesh: unknown keys {sorted(unknown)}")
        g = dict(DEPTH_MESH_GRADED)
        g.update(spec)
        if (
            not (0 < g["first_cell_cm"] <= g["near_ps_cm"])
            or g["ratio"] <= 1.0
            or g["junction_ps_cm"] <= 0
            or g["near_depth_cm"] <= 0
        ):
            raise ValueError(f"depth_mesh: invalid grading {g}")
        return g
    raise ValueError(f"depth_mesh must be 'graded', 'legacy' or a dict, got {spec!r}")


def junction_ps(spec=DEPTH_MESH_DEFAULT):
    """Spacing (cm) of the p+/n- junction line."""
    g = resolve_depth_mesh(spec)
    return LEGACY_JUNCTION_PS if g is None else g["junction_ps_cm"]


def contact_ps(spec=DEPTH_MESH_DEFAULT):
    """Spacing (cm) of the cathode contact line itself."""
    g = resolve_depth_mesh(spec)
    return LEGACY_CATHODE_PS if g is None else g["first_cell_cm"]


def contact_offsets(spec=DEPTH_MESH_DEFAULT, max_depth_cm=None):
    """Distances from the contact of the graded lines: ``[(d, ps), ...]``, d > 0
    increasing, ps the spacing of the cell beyond d. Empty for the legacy mesh.
    ``max_depth_cm`` (e.g. half the epi) truncates the list inside the layer."""
    g = resolve_depth_mesh(spec)
    if g is None:
        return []
    out = []
    d, h = 0.0, g["first_cell_cm"]
    while True:
        d += h
        h = min(h * g["ratio"], g["near_ps_cm"])
        if d >= g["near_depth_cm"] - 1e-12:
            break
        if max_depth_cm is not None and d >= max_depth_cm:
            break
        out.append((d, h))
    return out


def contact_lines(
    contact_pos_cm, direction, spec=DEPTH_MESH_DEFAULT, max_depth_cm=None
):
    """Absolute ``[(pos, ps)]`` of the contact line and its graded lines; ``direction``
    is -1 for lines BELOW the contact coordinate (a contact at the bottom of the layer,
    the cathode) and +1 for lines above it."""
    lines = [(contact_pos_cm, contact_ps(spec))]
    for d, ps in contact_offsets(spec, max_depth_cm):
        lines.append((contact_pos_cm + direction * d, ps))
    return lines

"""2D structures page: the fv2d engine (etna.core.fv2d) on its native meshes of the
mesa, guard-ring and axisymmetric 3D-electrode structures: bias, CCE vs position of a
normal alpha stripe, contact sharing, field / potential / doping maps.

All Streamlit calls live inside render(). The fv2d functions are called as module
attributes of etna.core.fv2d.structures / .workflows so tests can monkeypatch them. A
solve is 0.3-3 s (the first of a session compiles); results are cached per session.
The structures use the ALT_STACK (p+ anode, n- epi, ohmic cathode directly on the epi,
no n+ layer): they are NOT the sidebar device, so the page has its own inputs.
"""

from __future__ import annotations

import numpy as np
import plotly.graph_objects as go
import streamlit as st

from app.components import pubfig as P
from app.components.page_utils import about, show_error

STRUCTURES = {
    "mesa": "Mesa (pillar, trench side wall = zero-flux boundary)",
    "guard_ring": "Planar pad + guard ring (Cartesian half cell)",
    "electrode3d": "3D electrode, axisymmetric (r, z): n+ column, p+ wall",
}
MAX_BIAS_V = 200.0
# the impact-ionization sets of etna.core.physics_models are fitted from 1-2 MV/cm up
# (field_range_V_cm); fv2d has no avalanche term, so flag fields in that range
HIGH_FIELD_V_CM = 1.0e6


def _mods():
    from etna.core.fv2d import structures as fs
    from etna.core.fv2d import workflows as fw

    return fs, fw


def _geometry(kind):
    c = st.columns(4)
    g = {}
    if kind == "mesa":
        g["pillar_half_width_um"] = c[0].number_input(
            "Pillar half-width (µm)", 1.0, 500.0, 50.0, key="s2_pw"
        )
    elif kind == "guard_ring":
        g["sv_half_width_um"] = c[0].number_input(
            "Pad half-width (µm)", 1.0, 500.0, 50.0, key="s2_sv"
        )
        g["guard_ring_width_um"] = c[1].number_input(
            "Guard-ring width (µm)", 0.5, 100.0, 5.0, key="s2_gw"
        )
        g["guard_ring_gap_um"] = c[2].number_input(
            "Pad-to-guard gap (µm)", 0.2, 100.0, 3.0, key="s2_gap"
        )
        g["guard_ring_depth_um"] = c[3].number_input(
            "Guard-ring p+ depth (µm)", 0.05, 5.0, 1.0, key="s2_gd"
        )
        g["with_guard"] = st.checkbox("Guard ring present", True, key="s2_wg",
                                      help="Off: the same pad without the guard implant / contact.")  # fmt: skip
    else:
        g["outer_radius_um"] = c[0].number_input(
            "Cell outer radius (µm)", 10.0, 500.0, 50.0, key="s2_ro"
        )
        g["column_radius_um"] = c[1].number_input(
            "n+ column radius (µm)", 0.5, 100.0, 5.0, key="s2_rc"
        )
        g["top_gap_um"] = c[2].number_input(
            "Top gap around the column (µm)", 0.2, 100.0, 2.0, key="s2_tg"
        )
    return {k: (float(v) if not isinstance(v, bool) else v) for k, v in g.items()}


def solve_structure(p):
    """Build, bias and scan one structure (pure: no st.*)."""
    fs, fw = _mods()
    kw = dict(p["geometry"])
    kw.update(epi_thickness_cm=p["epi_um"] * 1e-4, T=p["T"])
    if p["kind"] != "electrode3d":
        kw["depth_mesh"] = p["depth_mesh"]
    builder = {"mesa": fs.mesa_structure, "guard_ring": fs.guard_ring_structure,
               "electrode3d": fs.electrode3d_structure}[p["kind"]]  # fmt: skip
    structure = builder(**kw)
    mesh = fs.build_mesh(structure)
    info = mesh.info
    solver = fw.prepared_solver(mesh, {"cathode": float(p["V"])})
    x_max = info.get("total_half_width_cm", info["half_width_cm"])
    xs = np.linspace(0.0, x_max, int(p["n"]))
    scans = {}
    if p["kind"] == "guard_ring" and info.get("guard_ring_contact"):
        ns = fw.cut_node_sets(mesh, info["sv_half_width_cm"] + 0.5 * info["gap_cm"])
        for c in ("anode", "guard_ring_anode"):
            scans[c] = fw.cce_lateral_scan_fv2d(
                solver, contact=c, node_sets=ns, x_positions_cm=xs
            )
    else:
        c = "cathode" if p["kind"] == "electrode3d" else "anode"
        scans[c] = fw.cce_lateral_scan_fv2d(solver, contact=c, x_positions_cm=xs)
    psi = np.asarray(solver.potential(), float)
    dpsi_dy, dpsi_dx = np.gradient(psi, mesh.ys, mesh.xs)
    E = np.hypot(dpsi_dx, dpsi_dy)  # V/cm
    return dict(
        xs_um=np.asarray(mesh.xs) * 1e4, ys_um=np.asarray(mesh.ys) * 1e4, psi=psi, E=E,
        net=np.asarray(mesh.net, float), scans=scans, x_scan_um=xs * 1e4,
        dark_check=getattr(solver, "dark_check_result", None), nx=mesh.nx, ny=mesh.ny,
        coordinate_system=mesh.coordinate_system,
        info={k: v for k, v in info.items() if np.isscalar(v) or v is None},
    )  # fmt: skip


def _inputs():
    kind = st.selectbox(
        "Structure", list(STRUCTURES), format_func=STRUCTURES.get, key="s2_kind"
    )
    geom = _geometry(kind)
    c = st.columns(4)
    epi = c[0].number_input("Epi thickness (µm)", 2.0, 200.0, 10.0, key="s2_epi")
    T = c[1].number_input("Temperature (K)", 200.0, 500.0, 300.0, key="s2_T")
    V = c[2].number_input(
        "Reverse bias (V, negative)", -MAX_BIAS_V, 0.0, -50.0 if kind == "electrode3d" else -10.0,
        key=f"s2_V_{kind}", help="App convention: reverse bias negative. fv2d applies +|V| to "
        f"the cathode with the anode(s) at 0 V. Capped at {MAX_BIAS_V:.0f} V.",
    )  # fmt: skip
    n = c[3].number_input("Scan positions", 3, 101, 21, key="s2_n",
                          help="Alpha-stripe positions from x = 0 to the cell edge.")  # fmt: skip
    dm = "graded"
    if kind != "electrode3d":
        dm = st.selectbox(
            "Depth mesh at the cathode", ["graded", "legacy"], key="s2_dm",
            help="graded (default since 2026-10-21): 0.1 nm cathode cell growing 1.25x; the "
            "ohmic cathode on the epi is a hole sink and the CCE converges at first order in "
            "that cell. legacy: the meshes before 2026-10-21 (100 nm cell; CCE +6e-3 at 10 V).",
        )  # fmt: skip
    return dict(kind=kind, geometry=geom, epi_um=float(epi), T=float(T), V=abs(float(V)),
                n=int(n), depth_mesh=dm)  # fmt: skip


def _scan_plot(r, p):
    xlab = "Radius r (µm)" if p["kind"] == "electrode3d" else "Lateral position x (µm)"
    fig = go.Figure()
    total = None
    for c, s in r["scans"].items():
        fig.add_trace(
            go.Scatter(
                x=r["x_scan_um"],
                y=s["cce_values"],
                mode="lines+markers",
                name=f"CCE at {c}",
            )
        )
        total = s["cce_values"] if total is None else total + s["cce_values"]
    if len(r["scans"]) > 1:
        fig.add_trace(
            go.Scatter(
                x=r["x_scan_um"],
                y=total,
                mode="lines",
                line_dash="dot",
                name="sum of contacts",
            )
        )
    fig.update_layout(xaxis_title=xlab, yaxis_title="CCE (collected / generated in the mesh)",
                      title="CCE vs alpha-stripe position")  # fmt: skip
    lo, hi = _y_range(r)
    if lo is not None:
        fig.update_yaxes(range=[lo, hi])
    return fig


def _y_range(r):
    """A y range of at least 0.01 so round-off (~1e-14 on a flat scan) is not drawn as
    structure; (None, None) when the data span more than that."""
    allv = np.concatenate(
        [np.asarray(s["cce_values"], float) for s in r["scans"].values()]
    )
    if allv.size and np.ptp(allv) < 1e-3:
        m = float(np.mean(allv))
        return m - 0.005, m + 0.005
    return None, None


def _map_plot(r, what):
    z, title, cs = {
        "|E| (kV/cm)": (r["E"] * 1e-3, "|E| (kV/cm)", "Viridis"),
        "Potential (V)": (r["psi"], "Potential (V)", "RdBu_r"),
        "log10 |net doping| (cm⁻³)": (
            np.log10(np.abs(r["net"]) + 1.0),
            "log10 |N_D − N_A| (cm⁻³)",
            "Plasma",
        ),
    }[what]
    kw = {}
    if what.startswith("|E|"):
        # the junction-edge peaks dominate the scale: clip it at the 99th percentile
        kw = dict(zmin=0.0, zmax=float(np.nanpercentile(z, 99.0)))
    fig = go.Figure(
        go.Heatmap(
            x=r["xs_um"],
            y=r["ys_um"],
            z=z,
            colorscale=cs,
            colorbar=dict(title=title),
            **kw,
        )
    )
    fig.update_layout(xaxis_title="Lateral x / radius r (µm)", yaxis_title="Depth y (µm)",
                      yaxis_autorange="reversed",
                      title=what + (f" (colour scale clipped at the 99th percentile, max "
                                    f"{np.nanmax(z):.0f} kV/cm)" if kw else ""))  # fmt: skip
    return fig


def _pub(r, p, what):
    prov = P.prov_for(solver="etna.core.fv2d (numba FV, devsim discretisation)")
    xlab = r"Radius r ($\mu$m)" if p["kind"] == "electrode3d" else r"Lateral x ($\mu$m)"
    notes = [f"{STRUCTURES[p['kind']]}; ALT_STACK (cathode on the epi, no n+); epi "
             f"{p['epi_um']:g} um, reverse bias {p['V']:g} V, {p['T']:g} K; mesh {r['nx']} x {r['ny']}."]  # fmt: skip

    def draw(fig):
        if what == "scan":
            ax = fig.add_subplot()
            for c, s in r["scans"].items():
                ax.plot(r["x_scan_um"], s["cce_values"], "o-", ms=2, label=c)
            ax.set_xlabel(xlab)
            ax.set_ylabel("CCE")
            lo, hi = _y_range(r)
            if lo is not None:
                ax.set_ylim(lo, hi)
                ax.ticklabel_format(axis="y", useOffset=False)
            ax.legend()
        else:
            ax = fig.add_subplot()
            im = ax.pcolormesh(
                r["xs_um"], r["ys_um"], r["E"] * 1e-3, shading="nearest", cmap="viridis",
                vmin=0.0, vmax=float(np.nanpercentile(r["E"] * 1e-3, 99.0)),
            )  # fmt: skip
            ax.invert_yaxis()
            ax.set_xlabel(xlab)
            ax.set_ylabel(r"Depth ($\mu$m)")
            ax.grid(False)
            fig.colorbar(im, ax=ax, label="|E| (kV/cm)", extend="max")

    if what != "scan":
        notes.append(
            f"Colour scale clipped at the 99th percentile of |E| (max "
            f"{np.nanmax(r['E']) * 1e-3:.0f} kV/cm)."
        )
    return lambda size: P.publication_figure(draw, prov, notes, size)


def render() -> None:
    st.title("2D Structures (fv2d)")
    st.caption(
        "Mesa, guard-ring and axisymmetric 3D-electrode structures on the general 2D "
        "engine: bias, CCE vs alpha position, contact sharing and field maps."
    )
    about(
        "- **Model:** `etna.core.fv2d`, a numba float64 box-method drift-diffusion solver "
        "with devsim's discretisation (Scharfetter-Gummel, SRH, ohmic contacts, static "
        "doping), on native tensor meshes; Cartesian or axisymmetric (r, z). Dark state by "
        "an adaptive Newton ramp, checked by a G = 0 Gummel pass (shown below).\n"
        "- **Cross-checked:** on devsim's mesh, parity with devsim is 1e-16 to 2.4e-12 in "
        "CCE and 1e-14 V in potential (mesa 50 V, guard ring 10 V incl. three-contact "
        "sharing, 3D electrode 50 V). Not validated against measurement.\n"
        "- **Device:** the ALT_STACK (p+ anode, graded n- epi, ohmic cathode directly on "
        "the epi, no n+ layer), NOT the sidebar device; epi thickness and temperature are "
        "set here. Incomplete ionization is not supported.\n"
        "- **CCE:** a normally incident Am-241 alpha stripe (Bragg profile, 2 µm wide) at "
        "each position; CCE = contact current / generation in the mesh. The mesa scan is "
        "flat by construction (zero-flux side wall, laterally uniform pillar).\n"
        "- **3D electrode:** the (r, z) cell with a p+ wall models annular-trench cells, "
        "NOT column-lattice 3D detectors (a square 50 µm cell depletes at 75 V vs 30 V for "
        "the equal-area cylinder; docs/USAGE.md 11.1).\n"
        "- **Mesh:** the graded cathode mesh is the default; the CCE converges at first "
        "order in the cathode cell (legacy mesh: +6e-3 at 10 V on the mesa)."
    )
    p = _inputs()
    tok = P.token(None, p)
    if st.button("Solve structure", key="s2_run", type="primary"):
        try:
            with st.spinner(
                "Building the mesh, ramping the bias and scanning (~1-5 s; the first run compiles)..."
            ):
                st.session_state["s2_result"] = P.cached_run(
                    "s2", tok, lambda: solve_structure(p)
                )
            st.session_state["s2_tok"] = tok
            st.session_state["s2_params"] = p
        except Exception as e:
            st.session_state.pop("s2_result", None)
            show_error(
                e,
                f"Try a lower bias (the ramp starts at 0 V; at most {MAX_BIAS_V:.0f} V) or the graded mesh.",
            )
    r = st.session_state.get("s2_result")
    if r is None:
        return
    if st.session_state.get("s2_tok") != tok:
        st.warning(P.STALE_TEXT)
    pu = st.session_state.get("s2_params", p)
    m = st.columns(4)
    m[0].metric("Mesh nodes", f"{r['nx']} x {r['ny']}")
    dc = r["dark_check"]
    if dc is not None:
        m[1].metric(
            "Dark check |du| (V)",
            f"{dc[1]:.1e}",
            help=f"{dc[0]} Gummel iterations at G = 0.",
        )
    first = next(iter(r["scans"].values()))
    m[2].metric("Max conservation residual", f"{np.nanmax(np.abs(first['residuals'])):.1e}",
                help="Sum of the balanced contact currents / I_gen (a real check with node-set cuts).")  # fmt: skip
    m[3].metric("In-mesh generation fraction", f"{first['in_mesh_generation_fraction']:.4f}",
                help="Fraction of the alpha's pairs generated inside the epi of the mesh.")  # fmt: skip
    st.caption(
        "PREDICTED (fv2d, cross-checked against devsim on its mesh; not validated against "
        "measurement)."
    )
    e_max = float(np.nanmax(r["E"]))
    if e_max > HIGH_FIELD_V_CM:
        st.warning(
            f"Peak |E| = {e_max / 1e6:.2f} MV/cm (at a junction edge or corner), inside the "
            "field range where the impact-ionization sets of the model library are measured "
            "(from 1-2 MV/cm); fv2d has no avalanche model. The CCE here assumes no multiplication; check the bias against "
            "the Breakdown & Avalanche page (1D parallel-plane, an upper bound)."
        )
    st.plotly_chart(_scan_plot(r, pu), key="s2_scan_plot")
    P.figure_block("s2_scan_fig", f"{pu['kind']}_cce_scan", _pub(r, pu, "scan"), st.session_state.get("s2_tok"))
    what = st.radio("Map", ["|E| (kV/cm)", "Potential (V)", "log10 |net doping| (cm⁻³)"],
                    horizontal=True, key="s2_map")  # fmt: skip
    st.plotly_chart(_map_plot(r, what), key="s2_map_plot")
    P.figure_block("s2_map_fig", f"{pu['kind']}_field_map", _pub(r, pu, "map"), st.session_state.get("s2_tok"))
    cols = {"position_um": r["x_scan_um"]}
    for c, s in r["scans"].items():
        cols[f"CCE_{c}"] = s["cce_values"]
        cols[f"residual_{c}"] = s["residuals"]
        cols[f"direct_mismatch_{c}"] = s["direct_mismatch"]
    hdr = [
        f"inputs: {pu}",
        "solver: etna.core.fv2d",
        "stack: ALT_STACK (cathode on the epi, no n+)",
    ]
    c = st.columns(2)
    c[0].download_button("Download scan CSV",
                         P.table_csv_bytes(f"fv2d {pu['kind']} CCE scan", cols, hdr),
                         f"fv2d_{pu['kind']}_scan.csv", "text/csv", key="s2_scan_csv", on_click="ignore")  # fmt: skip
    X, Y = np.meshgrid(r["xs_um"], r["ys_um"])
    c[1].download_button(
        "Download map CSV (one row per node)",
        P.table_csv_bytes(f"fv2d {pu['kind']} node maps",
                          {"x_um": X.ravel(), "y_depth_um": Y.ravel(), "potential_V": r["psi"].ravel(),
                           "E_abs_V_per_cm": r["E"].ravel(), "net_doping_cm3": r["net"].ravel()}, hdr),
        f"fv2d_{pu['kind']}_maps.csv", "text/csv", key="s2_map_csv", on_click="ignore",
    )  # fmt: skip

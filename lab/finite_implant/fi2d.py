"""Finite-p+-implant 2D cross-section (devsim, LAB ONLY): fringe capacitance at the p+ edge.

Geometry (x lateral, y depth, y = 0 the SiC top surface, cm internally):

    x = 0 .............. x_edge | x_edge+gap .................. x_max
    [ anode on p+ (y=0)       ] [ SiO2, y in [-t_ox, 0]          ]
                                [ field plate on the oxide top,  ]
                                [ x_edge+gap .. x_edge+L_fp, V=0 ]
    p+ (N_A, 0 <= y < y_j) only for x <= x_edge (abrupt lateral edge, planar)
    n- epi graded as DEFAULT_STACK; for x > x_edge the epi reaches the surface with
    N_D(y < y_j) = N_D_junction (the grading clamps at the junction value)
    n+ (1e19, 1 um) under the epi, cathode on its far face (y = total)
    fixed oxide charge Qf (cm^-2, positive) on the SiC surface nodes under the oxide

x = 0 and x = x_max are Neumann (mirror) walls: x = 0 is inside the implant
(the "plane" side), x_max is the far side of the termination.

Electrostatics: Poisson-only with two flat quasi-Fermi levels (reverse bias, negligible
current): n = n_i exp((psi - V_c)/V_t) (electrons in equilibrium with the cathode),
p = n_i exp((V_a - psi)/V_t) (holes with the anode). The SiC contacts are ohmic
(charge-neutral Dirichlet psi = V + V_t asinh(N/2n_i) via the etna/simple_physics
pattern), the field plate is an oxide contact (psi = V_a; phi_ms ignored / absorbed in Qf).

C per unit length (F/cm) = q d(sum NodeVolume * n)/dV_c (electrons enter only through
the cathode: this is cathode-to-(anode + plate), the measured C). Holes cross-check.
"""

import math

import devsim
import devsim.python_packages.simple_physics as sp
from devsim.python_packages.model_create import (
    CreateContactNodeModel,
    CreateContactNodeModelDerivative,
    CreateEdgeModel,
    CreateEdgeModelDerivatives,
    CreateNodeModel,
    CreateNodeModelDerivative,
    CreateSolution,
)
import numpy as np

from etna.core.device import DEFAULT_STACK as S
from etna.core.incomplete_ionization import ionized_acceptor_concentration
from etna.core.sic_material import SiC4H_Parameters, intrinsic_concentration

EPS0 = 8.854e-14
Q = 1.602e-19
KB_J = 1.3806503e-23
EPS_SIC = 9.7 * EPS0
EPS_OX = 3.9 * EPS0
CLAMP = 700


def _lines(mesh, d, pts):
    for pos, ps in sorted(pts):
        devsim.add_2d_mesh_line(mesh=mesh, dir=d, pos=pos, ps=ps)


def build(
    name="fi",
    x_edge_um=50.0,
    x_out_um=100.0,
    t_ox_um=1.0,
    fp_um=0.0,
    Qf=0.0,
    gap_um=0.05,
    edge_ps_um=0.02,
    far_ps_um=4.0,
    epi_um=10.0,
    T=300.0,
    null=0.0,
    surf_ps_um=0.01,
    NJ=None,
    NB=None,
):
    """Build the device; x_out_um = 0 -> laterally uniform full-anode reference (no oxide)."""
    p = SiC4H_Parameters()
    um = 1e-4
    yj = S["substrate_thickness_cm"]
    epi_end = yj + epi_um * um
    tot = epi_end + S["nplus_thickness_cm"]
    xe = x_edge_um * um
    xmax = xe + x_out_um * um
    fine = x_out_um > 0
    has_ox = fine and not null  # null: same mesh, p+ and anode over the whole width
    xa = xmax if null else xe  # anode / p+ extent
    t_ox = t_ox_um * um
    mesh = f"{name}_m"
    devsim.create_2d_mesh(mesh=mesh)

    # ---- x lines ----
    xo = xe + gap_um * um
    if fine:
        eps_ = edge_ps_um * um
        xs = [(0.0, far_ps_um * um), (xe, eps_), (xo, eps_), (xmax, far_ps_um * um)]
        if xe > 10 * um:
            xs.append((xe - 10 * um, 0.5 * um))
        xs.append((xo + 10 * um, 0.5 * um))
        if fp_um > 0:
            xfp = xo + fp_um * um
            xs += [(xfp, eps_)]
            if fp_um > 3:
                xs += [(xfp - 2 * um, 0.2 * um)]
            xs += [(xfp + 3 * um, 0.2 * um)]
        if x_out_um > 40:
            xs.append((xo + 40 * um, 2 * um))
        # de-duplicate / keep inside
        d = {}
        for pos, ps in xs:
            if 0 <= pos <= xmax + 1e-12:
                d[round(pos, 10)] = min(ps, d.get(round(pos, 10), 1.0))
        _lines(mesh, "x", list(d.items()))
    else:
        _lines(mesh, "x", [(0.0, xe), (xe, xe)])

    # ---- y lines (device2d pattern + extra epi lines) ----
    ys = [
        (0.0, surf_ps_um * um),
        (yj - 5e-6, 1e-6),
        (yj, 1e-7),
        (yj + 1e-4, 2e-6),
        (yj + 2e-4, 5e-6),
        (yj + 5e-4, 5e-6),
        (epi_end - 1e-4, 5e-6),
        (epi_end, 1e-7),
        (tot, 1e-5),
    ]
    AB = 1e-8  # thin buffer regions so devsim finds the contact nodes (device2d pattern)
    ys += [(-AB, AB), (tot + AB, AB)]
    if has_ox:
        ys.append((-t_ox, t_ox / 8))
        if fp_um > 0:
            ys.append((-t_ox - AB, AB))
    _lines(mesh, "y", ys)

    devsim.add_2d_region(
        mesh=mesh, material="SiC", region="sic", xl=0, xh=xmax, yl=0, yh=tot
    )
    devsim.add_2d_region(mesh=mesh, material="SiC", region="buf_a", xl=0, xh=xa, yl=-AB, yh=0)
    devsim.add_2d_region(
        mesh=mesh, material="SiC", region="buf_c", xl=0, xh=xmax, yl=tot, yh=tot + AB
    )
    if has_ox and fp_um > 0:
        devsim.add_2d_region(
            mesh=mesh, material="Ox", region="buf_fp", xl=xo, xh=xo + fp_um * um,
            yl=-t_ox - AB, yh=-t_ox,
        )
    if has_ox:
        devsim.add_2d_region(
            mesh=mesh, material="Ox", region="ox", xl=xo, xh=xmax, yl=-t_ox, yh=0
        )
        devsim.add_2d_interface(
            mesh=mesh,
            name="sic_ox",
            region0="sic",
            region1="ox",
            xl=xo,
            xh=xmax,
            yl=0,
            yh=0,
        )
    devsim.add_2d_contact(
        mesh=mesh,
        name="anode",
        material="metal",
        region="sic",
        xl=0,
        xh=xa,
        yl=0,
        yh=0,
        bloat=1e-10,
    )
    devsim.add_2d_contact(
        mesh=mesh,
        name="cathode",
        material="metal",
        region="sic",
        xl=0,
        xh=xmax,
        yl=tot,
        yh=tot,
        bloat=1e-10,
    )
    if has_ox and fp_um > 0:
        devsim.add_2d_contact(
            mesh=mesh,
            name="fplate",
            material="metal",
            region="ox",
            xl=xo,
            xh=xo + fp_um * um,
            yl=-t_ox,
            yh=-t_ox,
            bloat=1e-10,
        )
    devsim.finalize_mesh(mesh=mesh)
    devsim.create_device(mesh=mesh, device=name)

    # ---- SiC parameters & doping ----
    V_t = KB_J * T / Q
    n_i = intrinsic_concentration(T, p)[0]
    for k, v in (
        ("Permittivity", EPS_SIC),
        ("ElectronCharge", Q),
        ("n_i", n_i),
        ("V_t", V_t),
    ):
        devsim.set_parameter(device=name, region="sic", name=k, value=v)
    NA = ionized_acceptor_concentration(1e19, T)
    NJ = S["N_D_junction"] if NJ is None else NJ
    NB = S["N_D_bulk"] if NB is None else NB
    L = S["L_transition"]
    devsim.node_model(
        device=name,
        region="sic",
        name="Acceptors",
        equation=f"{NA}*step({yj}-y)*step({xe}-x)" if has_ox else f"{NA}*step({yj}-y)",
    )
    graded = f"({NB} + ({NJ} - {NB})*exp(-max(y - {yj}, 0)/{L}))"
    lat = f"max(step(y-{yj}), step(x-{xe}-1e-12))" if has_ox else f"step(y-{yj})"
    nplus = f" + {S['N_D_nplus']}*step(y-{epi_end})"
    devsim.node_model(
        device=name, region="sic", name="Donors", equation=f"{graded}*{lat}{nplus}"
    )
    # fixed oxide charge: sheet Qf on the y = 0 SiC nodes under the oxide (exact sheet
    # in the box discretisation: node volume = (h_y/2) * box width)
    CreateSolution(name, "sic", "Nfix")
    xn = np.array(devsim.get_node_model_values(device=name, region="sic", name="x"))
    yn = np.array(devsim.get_node_model_values(device=name, region="sic", name="y"))
    nfix = np.zeros_like(xn)
    if has_ox and Qf != 0.0:
        yl = np.unique(yn)
        h0 = yl[1] - yl[0]
        sel = (np.abs(yn) < 1e-12) & (xn >= xo - 1e-12)
        nfix[sel] = Qf / (h0 / 2)
    devsim.set_node_values(device=name, region="sic", name="Nfix", values=nfix.tolist())
    devsim.set_parameter(device=name, name="qfs", value=1.0)
    devsim.node_model(
        device=name, region="sic", name="NetDoping", equation="Donors-Acceptors+qfs*Nfix"
    )

    # ---- Poisson with two quasi-Fermi levels ----
    CreateSolution(name, "sic", "Potential")
    for c in ("anode", "cathode"):
        devsim.set_parameter(device=name, name=sp.GetContactBiasName(c), value=0.0)
    vc, va = sp.GetContactBiasName("cathode"), sp.GetContactBiasName("anode")
    en = f"n_i*exp(min(max((Potential-{vc})/V_t, -{CLAMP}), {CLAMP}))"
    hp = f"n_i*exp(min(max(({va}-Potential)/V_t, -{CLAMP}), {CLAMP}))"
    for nm, ex in (
        ("Electrons", en),
        ("Holes", hp),
        ("Rho", "kahan3(Holes, -Electrons, NetDoping)"),
        ("PotentialNodeCharge", "-ElectronCharge*Rho"),
    ):
        CreateNodeModel(name, "sic", nm, ex)
        CreateNodeModelDerivative(name, "sic", nm, ex, "Potential")
    for nm, ex in (
        ("ElectricField", "(Potential@n0-Potential@n1)*EdgeInverseLength"),
        ("PotentialEdgeFlux", "Permittivity*ElectricField"),
    ):
        CreateEdgeModel(name, "sic", nm, ex)
        CreateEdgeModelDerivatives(name, "sic", nm, ex, "Potential")
    devsim.equation(
        device=name,
        region="sic",
        name="PotentialEquation",
        variable_name="Potential",
        node_model="PotentialNodeCharge",
        edge_model="PotentialEdgeFlux",
        variable_update="log_damp",
    )
    # ohmic contacts: psi = V + V_t asinh(Net/2 n_i)  (charge neutral, QF = contact bias)
    for c in ("anode", "cathode"):
        b = sp.GetContactBiasName(c)
        eq = f"Potential - {b} - V_t*asinh(NetDoping/(2*n_i))"
        m = f"{c}nodemodel"
        CreateContactNodeModel(name, c, m, eq)
        CreateContactNodeModelDerivative(name, c, m, eq, "Potential")
        devsim.contact_equation(
            device=name,
            contact=c,
            name="PotentialEquation",
            node_model=m,
            edge_charge_model="PotentialEdgeFlux",
        )
    if has_ox:
        devsim.set_parameter(
            device=name, region="ox", name="Permittivity", value=EPS_OX
        )
        sp.CreateOxidePotentialOnly(name, "ox", "log_damp")
        sp.CreateSiliconOxideInterface(name, "sic_ox")
        if fp_um > 0:
            # plate tied to the anode bias (0 V): its own parameter, kept at va
            devsim.set_parameter(
                device=name, name=sp.GetContactBiasName("fplate"), value=0.0
            )
            sp.CreateOxideContact(name, "ox", "fplate")
    # initial guess: neutral potential
    net = np.array(
        devsim.get_node_model_values(device=name, region="sic", name="NetDoping")
    )
    psi0 = V_t * np.arcsinh(net / (2 * n_i))
    devsim.set_node_values(
        device=name, region="sic", name="Potential", values=psi0.tolist()
    )
    if has_ox:
        nox = len(devsim.get_node_model_values(device=name, region="ox", name="x"))
        devsim.set_node_values(
            device=name, region="ox", name="Potential", values=[0.0] * nox
        )
    info = dict(
        name=name,
        xe=xe,
        xmax=xmax,
        yj=yj,
        epi_end=epi_end,
        tot=tot,
        n_i=n_i,
        has_ox=has_ox,
        nodes=len(xn),
    )
    return info


def solve(rel=1e-11, it=60):
    devsim.solve(
        type="dc", absolute_error=1e30, relative_error=rel, maximum_iterations=it
    )


def set_vc(info, v):
    devsim.set_parameter(
        device=info["name"], name=sp.GetContactBiasName("cathode"), value=float(v)
    )


def _snap(info):
    regs = ["sic"] + (["ox"] if info["has_ox"] else [])
    return {r: devsim.get_node_model_values(device=info["name"], region=r, name="Potential")
            for r in regs}


def _restore(info, snap):
    for r, v in snap.items():
        devsim.set_node_values(device=info["name"], region=r, name="Potential", values=v)


def _step(info, v_from, v_to, depth=0):
    snap = _snap(info)
    set_vc(info, v_to)
    try:
        solve()
    except devsim.error:
        _restore(info, snap)
        if depth > 8:
            raise
        vm = 0.5 * (v_from + v_to)
        _step(info, v_from, vm, depth + 1)
        _step(info, vm, v_to, depth + 1)


def ramp(info, v_from, v_to, step=1.0):
    n = max(1, int(math.ceil(abs(v_to - v_from) / step)))
    vs = np.linspace(v_from, v_to, n + 1)
    for a, b in zip(vs[:-1], vs[1:]):
        _step(info, a, b)


def charges(info):
    nm = info["name"]
    vol = np.array(
        devsim.get_node_model_values(device=nm, region="sic", name="NodeVolume")
    )
    n = np.array(
        devsim.get_node_model_values(device=nm, region="sic", name="Electrons")
    )
    pp = np.array(devsim.get_node_model_values(device=nm, region="sic", name="Holes"))
    return Q * np.sum(vol * n), Q * np.sum(vol * pp)


def fields(info):
    nm = info["name"]
    g = lambda k: np.array(
        devsim.get_node_model_values(device=nm, region="sic", name=k)
    )
    return dict(
        x=g("x"),
        y=g("y"),
        n=g("Electrons"),
        p=g("Holes"),
        D=g("Donors"),
        psi=g("Potential"),
    )


def set_qfs(info, v):
    """Scale of the fixed oxide charge (continuation in Qf at fixed bias)."""
    devsim.set_parameter(device=info["name"], name="qfs", value=float(v))


def qf_continue(info, q_from, q_to, n=5, depth=0):
    """Step qfs q_from -> q_to at fixed bias, bisecting on failure."""
    qs = np.linspace(q_from, q_to, n + 1)
    for a, b in zip(qs[:-1], qs[1:]):
        snap = _snap(info)
        set_qfs(info, b)
        try:
            solve()
        except devsim.error:
            _restore(info, snap)
            set_qfs(info, a)
            if depth > 6:
                raise
            qf_continue(info, a, b, 2, depth + 1)

"""OPT-IN trap-model presets for the charged-trap / kinetic damage model (issue #5).

Nothing here is a default: ``DeviceConfig(radiation_traps=True)`` stays the Burin
Table I steady-state model (``radiation_traps.DEFAULT_TRAP_MODEL``). Each preset is a
``RadiationTrapModel`` with a ``label`` that says what is literature, what is a stated
assumption and what is CALIBRATED (and to which numbers). Evidence, fits and checks:
``lab/radiation_traps/notes/MECHANISMS.md`` and ``lab/radiation_traps/bench/mechanisms.py``.

Sources opened (full text unless marked):

- Z1/2 negative-U, two components: Capan et al., J. Appl. Phys. 124, 245701 (2018),
  arXiv:1812.06462, Table I (L-DLTS Arrhenius E_a and prefactor A per transition;
  [Z2]:[Z1] = 4.4 from long filling pulses). Each site is a two-level acceptor
  (0 -> - at the (-/0) level, - -> = at the deeper (=/-) level); the emission
  cross-section is converted from A to this code's prefactor convention
  (``arrhenius_sigma``), so the 300 K emission rates are Capan's (Z2(=/0) ~13 /s).
- Introduction rates Z1/2 5.0, EH4 2.4, EH6/7 1.6 cm^-1/neq: Burin et al.,
  arXiv:2407.16710 Table I (cited there to Hazdra et al. NIM B 327 (2014), not opened).
  The Z1/2 DLTS concentration is ASSUMED to count emitted electrons (single-level
  analysis of a two-electron emission): 2.5 V_C per neq with two electrons each
  (``z12_defects_per_neq``; 5.0 is the other end of that bracket).
- EH6/7 sigma_n: 1e-13 cm^2 at E = 1.60 eV, the value that reproduces the measured
  EH6/7 DLTS peak temperatures / LDLTS emission rates in this code's prefactor
  (Hazdra et al. SSP 242 (2016) 421 Fig. 1 and Table 1: Ec-1.64, 3e-13; Danno & Kimoto
  JAP 100 (2006) 113728; Alfieri & Kimoto APL 102 (2013) 152108; range 6e-14..2e-13).
  Burin/Gaggl's 9e-12 overshoots those emission rates by 40-160x and its cited source
  (Knezevic et al. Materials 16 (2023) 3347) contains no EH6/7 value; the audit-C-4
  1.4e-14 has no primary source and undershoots by 4-15x. The DEFAULT table keeps
  1.4e-14 (only K_tau and the steady depleted charge depend on it; see MECHANISMS.md).
- Missing removal: Radmanovac et al., arXiv:2510.11304 (252.7 MeV p, Table 2:
  g = 5.60 +- 0.81 onsemi 50 um 4e13, 4.60 +- 0.33 CNM W2 50 um, 4.48 +- 0.26 CNM W4
  100 um 4e14, per p/cm^2, C-V at 1 MHz). The table above gives about half of it in
  the depleted region (kinetic). The remainder is a component that stays charged in
  depletion AND counts once in the neutral bulk: at room temperature a slow deep
  acceptor (Lebedev et al. Materials 15 (2022) 8637: lower-half-gap acceptors; V_Si
  (0/-) near midgap, Bathen et al. npj Quantum Inf. 5 (2019), doi:10.1038/s41534-019-0227-y) and a direct donor loss are
  INDISTINGUISHABLE in C-V (same charge in both regions). The preset carries it as
  ``donor_removal`` (CALIBRATED); ``deep_acceptor=True`` puts the same rate on a
  midgap-like acceptor instead (differs only in CCE trapping and above ~500 K).
"""

from __future__ import annotations

import dataclasses

from etna.core.radiation_traps import (
    DEFAULT_TRAP_MODEL,
    DeepTrap,
    RadiationTrapModel,
    _vth,
)

# Capan et al. 2018, Table I: (E_a eV, A s^-1 K^-2) per transition
CAPAN_Z1 = {"-/0": (0.48, 7.4e6), "=/0": (0.59, 1.1e7)}
CAPAN_Z2 = {"-/0": (0.41, 2.15e6), "=/0": (0.67, 2.55e7)}
CAPAN_Z2_TO_Z1 = 4.4
# hole capture of Z1/2 (both levels): Burin Table I 3.5e-14 (unverified source)
Z12_SIGMA_P = 3.5e-14
EH67_SIGMA_N_EMISSION = 1e-13  # cm^2 at Ec-1.60 (see module docstring)

# CALIBRATED (lab/radiation_traps/bench/mechanisms.py, MECHANISMS.md): the slow
# removal component (cm^-1 per Si-neq, kappa(252.7 MeV) = 0.851 of this code) fitted
# to the three Radmanovac Table 2 slopes on top of literature_kinetic().
# Fit "d only" (3 slopes, chi2 1.13 for 2 dof): 2.87 +- 0.22 /neq (the data +- are the
# spread over depletion widths, not fit errors; kappa(252.7 MeV) = 0.851 of this code,
# 0.9 in the paper: a ~6 % systematic on the per-neq value, not on the per-p slope).
CALIBRATED_SLOW_REMOVAL = 2.87
CALIBRATED_SLOW_REMOVAL_SIGMA = 0.22


def arrhenius_sigma(A, T=300.0):
    """Emission cross-section (cm^2) in this code's convention, e_n = sigma v_th(T)
    N_C(T) exp(-E_a / kT), equal to an Arrhenius e_n = A T^2 exp(-E_a / kT)
    (v_th N_C ~ T^2: T only fixes the evaluation point)."""
    from etna.core.sic_material import SiC4H_Parameters, effective_dos

    NC, _ = effective_dos(float(T))
    return float(A) * float(T) ** 2 / (_vth(float(T), SiC4H_Parameters().m_e_dos) * NC)


def capan_z12_traps(defects_per_neq=2.5):
    """Z1 and Z2 (V_C at h and k) as two-level negative-U acceptors, Capan Table I
    energies and Arrhenius emission, split [Z2]:[Z1] = 4.4, total ``defects_per_neq``
    vacancies per neq (cm^-1)."""
    fz2 = CAPAN_Z2_TO_Z1 / (1.0 + CAPAN_Z2_TO_Z1)
    out = []
    for name, tab, frac in (("Z1", CAPAN_Z1, 1 - fz2), ("Z2", CAPAN_Z2, fz2)):
        E1, A1 = tab["-/0"]
        E2, A2 = tab["=/0"]  # rate-limiting first emission of the = state
        out.append(
            DeepTrap(
                name, "acceptor", defects_per_neq * frac, E1, arrhenius_sigma(A1),
                Z12_SIGMA_P, E2, arrhenius_sigma(A2), Z12_SIGMA_P,
            )  # fmt: skip
        )
    return tuple(out)


def literature_kinetic(z12_defects_per_neq=2.5, eh67_sigma_n=EH67_SIGMA_N_EMISSION,
                       lifetimes="occupancy", hold_s=60.0):  # fmt: skip
    """PREDICTED (no parameter fitted): negative-U Z1/Z2 (Capan), EH4 (Burin), EH6/7
    with the emission-consistent sigma_n, kinetic occupancy (0 V equilibrium start,
    1 V/s ramp, ``hold_s`` at bias for CCE), occupancy-weighted CCE lifetimes."""
    base = {t.name: t for t in DEFAULT_TRAP_MODEL.traps}
    traps = capan_z12_traps(z12_defects_per_neq) + (
        base["EH4"],
        dataclasses.replace(base["EH6/7"], sigma_n=float(eh67_sigma_n)),
    )
    return RadiationTrapModel(
        traps, occupancy="kinetic", hold_s=float(hold_s), lifetimes=lifetimes,
        label=(
            f"literature_kinetic: PREDICTED; Z1/Z2 negative-U (Capan 2018, "
            f"{z12_defects_per_neq:g} V_C/neq), EH4 + EH6/7 Burin Table I (EH6/7 "
            f"sigma_n {eh67_sigma_n:g}), kinetic occupancy"
        ),
    )  # fmt: skip


def deep_acceptor_trap(eta):
    """A midgap-like acceptor (Ec-1.9 eV: emission frozen at any lab temperature up to
    ~500 K) carrying a removal rate ``eta``. Its capture cross-sections are UNMEASURED
    (1e-14 both, assumed): they only enter the occupancy-weighted CCE lifetimes."""
    return DeepTrap("deep_acceptor", "acceptor", float(eta), 1.9, 1e-14, 1e-14)


def radmanovac_calibrated(deep_acceptor=False, slow_removal=None, **kw):
    """``literature_kinetic`` plus the slow removal component CALIBRATED to the three
    Radmanovac et al. Table 2 slopes (252.7 MeV p; ``CALIBRATED_SLOW_REMOVAL`` per neq
    unless ``slow_removal`` is given). ``deep_acceptor=False``: direct donor removal;
    True: the same rate on ``deep_acceptor_trap`` (identical room-temperature C-V)."""
    d = CALIBRATED_SLOW_REMOVAL if slow_removal is None else float(slow_removal)
    m = literature_kinetic(**kw)
    tag = (
        f"; + slow removal {d:.3g} /neq CALIBRATED to Radmanovac Table 2 "
        f"({'deep acceptor' if deep_acceptor else 'direct donor removal'})"
    )
    if deep_acceptor:
        return dataclasses.replace(
            m, traps=m.traps + (deep_acceptor_trap(d),), label=m.label + tag
        )
    return dataclasses.replace(m, donor_removal=d, label=m.label + tag)


PRESETS = {
    "burin_steady": lambda: DEFAULT_TRAP_MODEL,
    "burin_kinetic": lambda: dataclasses.replace(
        DEFAULT_TRAP_MODEL, occupancy="kinetic", label="burin_kinetic: Burin Table I"
    ),
    "literature_kinetic": literature_kinetic,
    "radmanovac_calibrated": radmanovac_calibrated,
    "radmanovac_calibrated_deep_acceptor": lambda: radmanovac_calibrated(True),
}


def trap_preset(name):
    """``RadiationTrapModel`` of preset ``name`` (``PRESETS``)."""
    if name not in PRESETS:
        raise ValueError(f"unknown trap preset {name!r} (have {sorted(PRESETS)})")
    return PRESETS[name]()

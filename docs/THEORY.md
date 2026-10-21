# ETNA theory and methods manual

This manual is for scientists who want to know **which equations ETNA actually solves**,
under which assumptions, with which constants, and how far each result can be trusted. Every
equation below was checked against the function that implements it (the file path and
function are given with it), and every constant against its source file. Where the code
departs from the textbook form, or a docstring disagrees with the code, a **Deviation**
or **Docstring note** says so.

What this manual is not:

- It is not the user guide. How to call the facades, presets, units of the API,
  solver selection, the app and the notebooks are in [docs/USAGE.md](USAGE.md), which this
  manual links to and does not repeat.
- It is not a validation report. The status labels below say what kind of evidence stands
  behind each ingredient; most outputs about the real detectors are predictions.

State of the code: `main` at the issue #9 follow-up merge (2026-10-21), package version
5.0.0. The defaults that changed on 2026-10-21 are the ones documented here:

- `DeviceConfig(radiation_traps="auto")` selects the calibrated kinetic trap preset;
- the NIEL convention is `"sic-damage"`;
- the low-field mobility set is `"c_axis"`.

Contents: [0 Conventions](#0-conventions) ·
[1 Device structure and doping](#1-device-structure-and-doping-profiles) ·
[2 Material parameters](#2-4h-sic-material-parameters) ·
[3 Electrostatics](#3-electrostatics-poisson-contacts-and-ionization) ·
[4 Drift-diffusion](#4-drift-diffusion-scharfetter-gummel-and-the-box-method) ·
[5 Mobility](#5-mobility) ·
[6 Recombination and dark current](#6-recombination-generation-and-the-dark-current) ·
[7 Impact ionisation](#7-impact-ionisation-and-breakdown-opt-in) ·
[8 Generation](#8-generation-by-radiation) ·
[9 Signal formation](#9-signal-formation-conduction-displacement-and-induced-current) ·
[10 CCE](#10-charge-collection-efficiency) ·
[11 FLASH](#11-flash-and-high-injection) ·
[12 Radiation damage](#12-radiation-damage) ·
[13 Temperature and self-heating](#13-temperature-and-self-heating) ·
[14 Microdosimetry](#14-microdosimetry) ·
[15 Monte Carlo deposition](#15-monte-carlo-energy-deposition-and-straggling) ·
[16 Numerics](#16-numerics) ·
[17 Uncertainty](#17-uncertainty) ·
[References](#references)

**Worked examples.** Each chapter ends with a short example that was run on this code; all of them
are in `lab/theory_manual/worked_examples.py` (one function per chapter, about 10 s in total after
numba compilation; the first run compiles the track engine for about 30 s):

```bash
uv run python lab/theory_manual/worked_examples.py              # all chapters
uv run python lab/theory_manual/worked_examples.py niel damage  # a subset
```

The output of the run quoted here is in `lab/theory_manual/notes/LOG.md`.
`tests/test_theory_manual.py` checks that every `etna/...` path and every function named in this
manual exists, so the manual cannot silently drift from the code.

## 0. Conventions

### 0.1 Status labels

Each ingredient carries one of these labels. They follow the "What to trust" table of the
README and [USAGE §2](USAGE.md#2-what-to-trust), and a label is never stronger than the one
given there.

| label | meaning | examples |
| --- | --- | --- |
| **verified** | code-to-code agreement with devsim, or with an analytic limit, inside a stated box | numba engines vs devsim; SG flux vs exact drift; Hecht closed form |
| **cross-checked** | compared with a literature result on *another* device or with an independent dataset; says the scale is right, not that this detector is modelled | C-V scale vs Bruzzi & Verroi 2023 (a Schottky diode on a different wafer); breakdown vs Hatakeyama's mesas |
| **calibrated** | a parameter fitted to data, so agreement with that data is fit quality, not validation | the doping profile (one C-V curve); the leakage scale $N_t$ (one upper bound); the direct donor removal (three 252.7 MeV C-V slopes of *other* devices) |
| **literature value** | an input constant read from a cited source (primary or a review) | mobility, Auger, impact ionisation, trap levels, NIST stopping powers |
| **assumption** | a modelling choice without a source, stated as such | the legacy Al ionized fraction; $f_{\rm depl}$ of the dark-current model; $N_t(T)\propto n_i(T)$ |
| **predicted** | a model output not checked against a measurement of this device | CCE vs bias, CCE vs fluence, temperature shapes, FLASH sensitivity, PHD |

### 0.2 Units and signs

The code works in CGS-like device units: lengths in cm, potentials in V, time in s,
densities in cm$^{-3}$, mobilities in cm$^2$/(V s), energies of levels and gaps in eV,
temperatures in K, currents in A/cm$^2$ (1D) or A/cm (2D, per unit depth). User-facing
inputs use µm and keV/µm (LET); NIEL and stopping powers are mass quantities in
MeV cm$^2$/g. The formulas below use these units throughout.

The 1D diode has the p$^+$ anode at $x=0$ (the irradiated surface), the junction at
$x_j$ = the p$^+$ thickness, the n$^-$ epi to $x_j+t_{\rm epi}$, and the n$^+$ layer behind
it with the cathode on its far end. Reverse bias is a negative number in the API; internally
the cathode is driven to $+|V|$ ([USAGE §5.4](USAGE.md#5-the-device-stack-presets-units-and-sign-conventions)).
The potential $\psi$ is referred to the intrinsic level (equilibrium Fermi level $=0$).
Under reverse bias electrons drift to $+x$ (cathode) and holes to $-x$ (anode).

### 0.3 Physical constants as coded

| symbol | value in the code | where |
| --- | --- | --- |
| $q$ | $1.602\times10^{-19}$ C | `etna/core/device.py` `Q`, `etna/core/analytical.py`, `etna/core/batched/material.py` |
| $\varepsilon_0$ | $8.854\times10^{-14}$ F/cm | `etna/core/device.py` `EPS_0` |
| $k_B$ | $1.3806503\times10^{-23}$ J/K (solvers' $V_t$); $8.617\times10^{-5}$ eV/K (`compute_ni`, `analytical.py`, `dark_current.py`) | `etna/core/device.py` `K_B_J` |
| $V_t = k_BT/q$ | 0.025855 V at 300 K | `etna/core/device.py` |

**Deviation.** The two Boltzmann constants give $V_t = 0.025855$ V in the solvers but
$k_BT = 0.025851$ eV in `compute_ni` and the analytic formulas, a $1.5\times10^{-4}$ relative
inconsistency; it is far below every other uncertainty and is kept so that pinned numbers
stay bit-identical.

### 0.4 Symbol table

| symbol | meaning | unit |
| --- | --- | --- |
| $\psi$ | electrostatic potential (intrinsic-level reference) | V |
| $\phi_n,\ \phi_p$ | electron / hole quasi-Fermi potentials | V |
| $u = \psi/V_t$ | scaled potential of the numba engines | 1 |
| $n,\ p$ | electron / hole densities | cm$^{-3}$ |
| $n_i$ | intrinsic density (anchored, section 2.2) | cm$^{-3}$ |
| $N_C,\ N_V$ | effective densities of states | cm$^{-3}$ |
| $N_D,\ N_A$ | donor / acceptor densities; $N_D^+,\ N_A^-$ ionized | cm$^{-3}$ |
| $N_{\rm eff}$ | effective (net space-charge) doping | cm$^{-3}$ |
| $\varepsilon = \varepsilon_r\varepsilon_0$ | permittivity, $\varepsilon_r = 9.7$ | F/cm |
| $E_g$ | band gap | eV |
| $\mu_n,\ \mu_p$ | mobilities | cm$^2$/(V s) |
| $v_{\rm sat}$ | saturation velocity | cm/s |
| $D = \mu V_t$ | diffusivity | cm$^2$/s |
| $\mathbf J_n,\ \mathbf J_p$ | current densities | A/cm$^2$ |
| $U$ | net recombination rate | cm$^{-3}$ s$^{-1}$ |
| $G$ | generation rate (radiation or effective) | cm$^{-3}$ s$^{-1}$ |
| $\tau_n,\ \tau_p$ | SRH lifetimes | s |
| $C_n,\ C_p$ | Auger coefficients | cm$^6$/s |
| $\alpha_n,\ \alpha_p$ | impact-ionisation coefficients | cm$^{-1}$ |
| $B(x) = x/(e^x-1)$ | Bernoulli function | 1 |
| $h,\ V_k,\ S_{kl}$ | edge length, control volume, face area of the box method | cm, cm$^3$ (1D: cm), cm$^2$ |
| $W$ | depletion width; $W_C = \varepsilon A/C$ its capacitance definition | cm |
| $C$ | capacitance | F |
| $A$ | device area (4 mm$^2$ = 0.04 cm$^2$ dosimeter) | cm$^2$ |
| $t_{\rm epi},\ x_j$ | epi thickness, junction depth (= p$^+$ thickness) | cm |
| $E_{\rm pair}$ | e-h pair creation energy, 7.8 eV | eV |
| $S/\rho,\ L$ | mass and linear stopping power (LET) | MeV cm$^2$/g, keV/µm |
| $R$ | CSDA range | cm |
| $\Phi,\ \Phi_{\rm eq}$ | particle fluence, equivalent (neq) fluence | cm$^{-2}$ |
| $\kappa(E)$ | NIEL hardness factor (chapter 12) | 1 |
| $\kappa_{\rm t}(E)$ | SiC-to-tissue LET ratio (chapter 14) | 1 |
| $\eta_k$ | defect introduction rate per neq | cm$^{-1}$ |
| $\sigma_n,\ \sigma_p$ | capture cross-sections | cm$^2$ |
| $v_{\rm th} = \sqrt{3k_BT/m^*}$ | thermal velocity | cm/s |
| $e_n,\ e_p$ | emission rates | s$^{-1}$ |
| $f$ | electron occupancy of a trap level | 1 |
| $\Gamma$ | Hurkx field-enhancement factor | 1 |
| $M$ | avalanche multiplication | 1 |
| $y,\ \bar l$ | lineal energy, mean chord length | keV/µm, µm |
| $y_F,\ y_D$ | frequency- and dose-mean lineal energy | keV/µm |
| CCE | charge collection efficiency | 1 |

## 1. Device structure and doping profiles

### 1.1 The stacks

The stack records are dictionaries in `etna/core/device.py`; `DeviceConfig`
(`etna/api/device.py`) converts its µm fields to them (`device_kwargs`). The default since
2026-10-18 is `DEFAULT_STACK` (alias `MEASURED_CV_STACK`):

| layer | default | note |
| --- | --- | --- |
| p$^+$ anode | 0.3 µm, $N_A = 10^{19}$ cm$^{-3}$ (Al) | `DeviceConfig.substrate_thickness_um` is the **p$^+$ anode**, despite its name (USAGE §5.1) |
| n$^-$ epi | 10 µm, graded profile (1.3) | `N_D_junction` 8.185880e15, `N_D_bulk` 3.382506e13, `L_transition` 0.4638630 µm |
| graded n$^-$/n$^+$ transition | last 1.136634 µm of the epi, to `buffer_end_doping` 4.516421e16 | log-linear in $N_D$ |
| n$^+$ | 1 µm, $10^{19}$ cm$^{-3}$ | truncated substrate, cathode on its far end |
| area | 0.04 cm$^2$ (4 mm$^2$ dosimeter) | `AREA_DOSIMETER_CM2` |

Other records: `NPLUS_CALIBRATION` (the previous default, abrupt n$^+$ step),
`LEGACY_STACK` (p$^+$ 1 µm, no n$^+$, cathode on the epi), and the opt-in
`MEASURED_CV_TABULATED_STACK` (1.4). The microdosimeter preset
(`MICRODOSIMETER_PRESET`, `etna/api/device.py`) is a 10 µm epi with uniform
$N_D = 10^{14}$, p$^+$ 0.3 µm, n$^+$ 1 µm, area 0.25 cm$^2$, no transition layer
(geometry from Petringa et al. 2025 [P25a, P25b]).

The whole stack is **one** SiC region: p$^+$, epi and n$^+$ differ only in their doping,
so there is no heterointerface and no interface charge. There are no oxides, no surface
states and no edge termination on the default path; the 2D structures with finite p$^+$,
mesa or guard ring are opt-in (`etna/core/alternative_structures.py`, USAGE §7).

**1D mesh** (`build_diode_mesh`, `etna/core/device.py`): mesh lines with spacing 100 nm at
the anode, 10 nm at $x_j-50$ nm, 1 nm at the junction $x_j$, 50 nm at $x_j+2$ µm and
$x_j+5$ µm, 1 nm at the epi end when an n$^+$ layer follows (100 nm without), and 100 nm in
the n$^+$; devsim grades between these lines. The default device has 684 nodes, 255 of them in
the 1.14 µm transition (spacing 1-12 nm). The batched engine reads devsim's node positions
(`etna/core/batched/mesh.py`), so both solvers discretise on the same nodes.

### 1.2 Acceptors

$$
N_A^{\rm node}(x) = f_A\,N_A\;\mathrm{step}(x_j-x),\qquad \mathrm{step}(0)=1 ,
$$

(`set_graded_doping_profile`, `set_doping_profile`, `set_tabulated_doping_profile`). The
devsim node model `Acceptors` holds the **ionized** density: $f_A = 0.132222$ at
$10^{19}$, 300 K with the default ionization model (3.4). Because `step(0)` = 1 the junction
node carries both dopants. The total $N_A$ (for the mobility) is recovered as
$N_A^{\rm node}/f_A$ (`acceptor_total_scale`).

### 1.3 The graded epi and the n$^-$/n$^+$ transition

With $d = x-x_j$ the depth into the epi, the base profile is (`graded_base_expr`; numpy twin
`donor_profile`)

$$
N_{\rm base}(x) = N_b + (N_j - N_b)\,\exp\!\left(-\frac{\max(d,0)}{L}\right),
$$

with $N_j$ = `N_D_junction`, $N_b$ = `N_D_bulk`, $L$ = `L_transition`. "uniform" is
$N_j=N_b=N_D$. The transition layer of thickness $L_b$ ending at $x_{\rm end}=x_j+t_{\rm epi}$,
$x_0 = x_{\rm end}-L_b$, is (`_buffer_term`, `buffer_spec`)

$$
t = \min\!\big(\max((x-x_0)/L_b,\,0),\,1\big),\qquad
N_{\rm gr} = \exp\!\big((1-t)\ln N_{\rm base} + t\ln N_{\rm end}\big),
$$

$$
N_D(x) = N_{\rm base}\,\mathrm{step}(d) + \mathrm{step}(x-x_0)\,\big(1-\mathrm{step}(x-x_{\rm end})\big)\big(\max(N_{\rm base},N_{\rm gr})-N_{\rm base}\big) + N_{D,n^+}\,\mathrm{step}(x-x_{\rm end}) ,
$$

so the transition is log-linear from the local epi value to $N_{\rm end}$ and never lowers
the epi value. **Deviation** (minor, stated in `lab/buffer_profile/notes/LOG.md`): at the node
$x = x_{\rm end}$ the donors are $N_{\rm base}+N_{D,n^+}$, not $N_{\rm end}+N_{D,n^+}$, because
the transition switches off where the n$^+$ step switches on; the difference is below
$10^{-2}$ of $N_{D,n^+}$.

### 1.4 The tabulated profile (opt-in)

`DeviceConfig(donor_table=...)` (preset `MEASURED_CV_TABULATED_PRESET`, table
`MEASURED_CV_TABLE`: 19 knots from 0 to 10 µm) interpolates linearly in $\log_{10}N_D$:

$$
\log_{10}N_D(x) = \mathrm{interp}\big(\mathrm{clip}(x-x_j,\,0,\,d_{\rm last});\ d_k,\ \log_{10}N_k\big),\quad x\ge x_j ,
$$

plus the n$^+$ step (`donor_table_spec`, `tabulated_donor_profile`). The values are written
into devsim as a node *solution*, so devsim and both numba engines see bit-identical floats.
Details, gates and parity are in USAGE §5 and `lab/tabulated_profile/notes/LOG.md`.

### 1.5 The n$^+$ substrate

The real substrate (hundreds of µm) is replaced by 1 µm of $10^{19}$ n$^+$ with the cathode on
its far end; $10^{18}$ and $10^{19}$ give the same $W$, CCE and dark current to $<10^{-3}$
(`etna/core/device.py` comment at `DEFAULT_N_D_NPLUS`). Radiation generation beyond the epi
end is dropped by default (`confine_generation_to_epi`, `nplus_generation=False`), a
modelling choice for the alpha measurements whose tracks do not reach the substrate.
`nplus_thickness_um=None` restores the legacy cathode on the epi; that model has a
convergence wall where the depletion reaches the contact, a contact artefact, not
punch-through (section 16.5).

**Docstring notes.** `create_sic_device` says "the n$^+$ inherits the device's scalar
mobility": since 2026-10-18 the mobility is local per node (5.1); only the lifetimes are still
scalars. `DeviceConfig` describes the n$^+$ as "fully ionized": true only for the default
ionization model; with `incomplete_ionization=True` the neutral $10^{19}$ n$^+$ is 51 % ionized.

### 1.6 Where the default profile comes from

**Status: calibrated.** The `DEFAULT_STACK` doping is a fit of the model's
$\varepsilon A/C$ to 230 points of the digitised C-V curve of Fig. 4b of the "Photons" paper
[PHOT] (0 to $-60$ V, $A$ = 4 mm$^2$; `lab/cv_gap/notes/LOG.md` §1-6): $\chi^2$ 198/224,
plateau rms 0.013 µm, knee rms 0.24 µm. Agreement with that curve is fit quality, not
validation; there is no SIMS, and the fitted bulk $3.4\times10^{13}$ cm$^{-3}$ is below the
documented 0.5-1$\times10^{14}$. The knee gate "within the digitisation uncertainty" fails for
the parametric form; the tabulated profile meets it ($\chi^2$ 30.8/224, knee rms 0.015 µm).
Gates: `tests/test_measured_cv_preset.py`.

**Worked example** (`worked_examples.py device`, section 5 of its output): `run_cv(DeviceConfig())`
gives $\varepsilon A/C$ = 1.6547 / 2.5402 / 4.1724 / 6.9390 / 9.2520 / 9.5185 µm at
0 / -0.5 / -1 / -2 / -5 / -10 V, reach-through 3.1 V, and the analytic full-depletion voltage
(moment integral of 3.6 up to the transition start at 8.863 µm) is 2.736 V.

## 2. 4H-SiC material parameters

### 2.1 Constants

`SiC4H_Parameters` (`etna/core/sic_material.py`) is the single source; other modules read it.

| quantity | value | source in the code | used by the solvers |
| --- | --- | --- | --- |
| $E_g(300)$ | 3.26 eV | Ioffe NSM [IOFFE] | via `bandgap` |
| Varshni $E_g(0),\ \alpha,\ \beta$ | 3.2965625 eV, $6.5\times10^{-4}$ eV/K, 1300 K | $\alpha,\beta$ Ioffe; $E_g(0)$ chosen so that $E_g(300)=3.26$ | yes |
| $\varepsilon_r$ | 9.7 (T-independent) | Ioffe | yes |
| $m_{e,\rm dos},\ m_{h,\rm dos}$ | $0.77\,m_0$, $1.0\,m_0$ | Ioffe ($m_e$ includes the 3 valleys) | $N_C, N_V$, $v_{\rm th}$ |
| $n_i(300)$ | $5.0\times10^{-9}$ cm$^{-3}$ | Ioffe / TU Wien [AYA] | the anchor of $n_i(T)$ |
| $\tau_n,\ \tau_p$ (300 K) | $10^{-9}$ s, $6\times10^{-7}$ s | Ioffe | yes |
| $\alpha_\tau$ | 1.72 | Ayalew (Z1/2) | $\tau(T)$ |
| $C_n,\ C_p$ | $5\times10^{-31}$, $2\times10^{-31}$ cm$^6$/s | Galeckas et al. 1997 [GAL97] | when Auger is on |
| $E_t$ (Z1/2), $m_t$ | $E_c-0.67$ eV, $0.25\,m_0$ | Burin Table I [BUR24t]; $m_t$ assumed | dark-current model |
| $N_t$ | $1.355575\times10^{12}$ cm$^{-3}$ s$^{-1}$ | **calibrated** (6.3) | dark-current model |
| $\rho$ | 3.21 g/cm$^3$ | Ioffe | generation, stopping, thermal |
| $E_{\rm pair}$ | 7.8 eV | "Lebedev; ~7.7-7.8 eV" (no full reference in the code) | generation |
| $B_{\rm rad}$ | $1.5\times10^{-12}$ cm$^3$/s | estimate | **not used** |
| $E_A,\ g_A$; $E_{D,h},\ E_{D,c}$ | 0.220 eV, 4; 0.050, 0.092 eV | TU Wien | **not read** by any solver (3.4) |

**Docstring note.** The comment at `E_D_hex` says donors are "taken fully ionized, the
incomplete-ionization model covers the Al acceptor only". This is stale: the opt-in model of
3.4 ionizes N donors too, with Kagamihara's energies (70.9 / 123.7 meV), not 50 / 92 meV.

### 2.2 Temperature laws

$$
E_g(T) = 3.2965625 - \frac{6.5\times10^{-4}\,T^2}{T+1300}\ \text{eV},
$$

$$
N_{C,V}(T) = 2\left(\frac{2\pi m_{e,h}k_BT}{h^2}\right)^{3/2},\qquad
n_{i,\rm raw}(T) = \sqrt{N_CN_V}\,e^{-E_g/2k_BT},
$$

$$
n_i(T) = n_{i}(300)\,\frac{n_{i,\rm raw}(T)}{n_{i,\rm raw}(300)},
$$

in `bandgap`, `compute_ni` and `intrinsic_concentration` (`etna/core/sic_material.py`);
$\tau(T) = \tau_{300}(T/300)^{1.72}$ (`srh_lifetime`), and the mobility laws of chapter 5.

**Deviation.** The solvers use the anchored $n_i$, $5.0\times10^{-9}$ at 300 K, whereas the
first-principles value from the same $N_C$, $N_V$ and $E_g$ is $8.52\times10^{-9}$ (factor
0.587, i.e. an effective gap 3.2876 eV at 300 K). `effective_dos` returns the *unanchored*
$N_C,N_V$; the opt-in incomplete ionization and the trap occupancy (chapter 12) use them, so
$n_i^2 \ne N_CN_Ve^{-E_g/k_BT}$ there by a factor 0.344. The trap module ties $p_1 = n_i^2/n_1$
to the solver's $n_i$ to keep detailed balance (12.5).

Not modelled: $\varepsilon_r(T)$, band-gap narrowing (tabulated, not coupled: 5.4), electron
affinity and Schottky barriers (all contacts are ohmic, 3.2).

**Status:** inputs are **literature values**; $E_g(0)$ and the $n_i$ anchor are ETNA choices;
the 300 K pins are tested (`tests/test_material.py`); every temperature *shape* is
**predicted** (an assumption-driven extrapolation).

**Worked example** (`worked_examples.py materials` and `device`): $E_g$ = 3.260000 / 3.248305 /
3.235386 eV at 300 / 350 / 400 K; $N_C = 1.6951\times10^{19}$, $N_V = 2.5088\times10^{19}$
cm$^{-3}$ (tabulated 1.69e19 / 2.49e19, pinned within 1 %); $n_{i,\rm raw}(300) = 8.5207\times10^{-9}$;
solver $n_i$ = 9.98e-15 / 5.0e-9 / 6.25e-5 / 7.72e-2 cm$^{-3}$ at 250 / 300 / 350 / 400 K.

## 3. Electrostatics: Poisson, contacts and ionization

### 3.1 Poisson equation

$$
\nabla\cdot(\varepsilon\nabla\psi) = -q\,\big(p - n + N_D^+ - N_A^- + \rho_{\rm trap}/q\big),
$$

discretised with the box method of chapter 4 (edge flux $\varepsilon(\psi_k-\psi_l)/h$).

- **Equilibrium (Poisson only)**, `_create_sic_potential_only` (`etna/core/poisson.py`):
  Boltzmann carriers $n = n_i\exp(\psi/V_t)$, $p = n_i\exp(-\psi/V_t)$ with the exponent
  clamped at $\pm700$; charge $p-n+$`NetDoping` summed with `kahan3`; update `log_damp`.
  Poisson-only under bias is not used (audit A-06): every biased state goes through DD.
- **Drift-diffusion**, `setup_sic_drift_diffusion` (`etna/core/drift_diffusion.py`): $n$, $p$
  become solution variables and the charge is `kahan3(Holes, -Electrons, NetDoping)`.
- **Batched engine**, `etna/core/batched/poisson.py`: the quasi-Fermi form
  $n = n_ie^{(\psi-\phi_n)/V_t}$, $p = n_ie^{(\phi_p-\psi)/V_t}$ with **flat** quasi-Fermi
  levels $\phi_n = V_{\rm cathode}$, $\phi_p = V_{\rm anode}$ for C-V, $W$ and fields. In
  4H-SiC ($n_i\approx5\times10^{-9}$) the reverse current is negligible, so this *is* the DD
  electrostatics (verified below). Scaled residual, $u=\psi/V_t$, $N_0 = \max|N_{\rm net}|$,
  $\lambda^2 = \varepsilon V_t/(qN_0)$:

$$
F_i = \lambda^2\left[\frac{u_{i+1}-u_i}{h_i}-\frac{u_i-u_{i-1}}{h_{i-1}}\right] + v_i\,(\hat p-\hat n+\hat D)_i ,
$$

  Newton with a Thomas solve, the update clamped at 20 $V_t$ (2 $V_t$ with traps) and the
  exponents at 60.

$\rho_{\rm trap}$ is non-zero only with the charged radiation traps (12.5).

### 3.2 Ohmic contacts

Contacts are charge-neutral ohmic Dirichlet conditions (devsim
`simple_physics.CreateSiliconPotentialOnlyContact` / `CreateSiliconDriftDiffusionAtContact`,
restated byte for byte in `contact_values`, `etna/core/batched/poisson.py`). With the net
doping $N$ at the contact node,

$$
n_c = 10^{-10} + \tfrac12\left|N+\sqrt{N^2+4n_i^2}\right|,\qquad
p_c = 10^{-10} + \tfrac12\left|-N+\sqrt{N^2+4n_i^2}\right|,
$$

$$
\psi_c = V_{\rm bias} + V_t\ln(n_c/n_i)\ (N>0),\qquad \psi_c = V_{\rm bias} - V_t\ln(p_c/n_i)\ (N\le0),
$$

and the DD contact rows fix $n=n_c$, $p=n_i^2/n_c$ (n-type) or the mirror (p-type), i.e.
infinite surface recombination velocity. With incomplete ionization $N$ is the neutral-state
ionized net doping `IonNetDopingEq`.

### 3.3 Built-in potential

The solver's built-in voltage is the contact difference. With the n$^+$ layer it contains the
p$^+$/n$^-$ junction *and* the n$^-$/n$^+$ high-low step:

$$
V_{bi}^{\rm solver} = V_t\ln\frac{N_A^-\,N_{D,n^+}}{n_i^2} = 3.1983\ \text{V}\quad(\psi_{\rm anode} = -1.5730,\ \psi_{\rm cathode} = +1.6253\ \text{V}).
$$

The analytic formula `built_in_potential` (`etna/core/analytical.py`) is the junction alone,
$V_{bi} = k_BT\ln(N_A^-N_{D}/n_i^2)$ with the ionized $N_A^-$: 3.0141 V with
$N_D = N_j$. **Docstring note:** it says "~15 % ionization, $N_A^-\approx1.5\times10^{18}$";
the code uses 13.2 % ($1.322\times10^{18}$).

### 3.4 Incomplete ionization

**Default (`incomplete_ionization=False`): a fixed ionized fraction of the Al acceptor**
(`ionized_acceptor_fraction`, `etna/core/incomplete_ionization.py`):

$$
E_{A,\rm eff} = \max(0.220 - 3.4\times10^{-8}N_A^{1/3},\,0),\quad
f_{\rm G} = \frac{1}{1+g_A e^{E_{A,\rm eff}/k_BT}},\quad g_A = 4,
$$

$$
f_{\rm emp} = \mathrm{clip}\big(0.05+0.10(\log_{10}N_A-18),\,0.01,\,0.95\big)\,\max\!\big(1+0.5\,\tfrac{T-300}{300},\,0.5\big),
$$

$$
w = \frac{1}{1+e^{-(\log_{10}N_A-18)/0.5}},\qquad f_A = \mathrm{clip}\big((1-w)f_{\rm G}+wf_{\rm emp},\,0,\,1\big).
$$

At $10^{19}$, 300 K: $f_A = 0.132222$ ($N_A^- = 1.3222\times10^{18}$); 0.1546 at 400 K. The
fraction is a scalar per device, independent of bias and position; donors are fully ionized.
**Status: assumption.** The coefficient $3.4\times10^{-8}$ and the empirical branch are
uncited; the module cites Ayalew §3.6, Pensl & Choyke 1993 and Kimoto & Cooper, but no value
is traceable to them. The low-doping branch rises with $N_A$ (audit-v7 A-02), so
`DeviceConfig` rejects $N_A<10^{17}$.

**Opt-in (`incomplete_ionization=True`): carrier-dependent ionization in Poisson**, the
steady Boltzmann occupancy of the review [BUR24r, Eqs. 82, 85, 88]:

$$
N_D^+ = \sum_l \frac{w_lN_D}{1+g_D\,n/n_{1,l}},\quad n_{1,l} = N_Ce^{-\Delta E_l/V_t};\qquad
N_A^- = \sum_l\frac{w_lN_A}{1+g_A\,p/p_{1,l}},\quad p_{1,l} = N_Ve^{-\Delta E_l/V_t},
$$

$$
\Delta E(N) = \max(\Delta E_0-\alpha N^{1/3},\,0).
$$

| dopant | levels $\Delta E_0$, $\alpha$, share | $g$ | source (opened in the review) |
| --- | --- | --- | --- |
| N | hexagonal 70.9 meV, $3.38\times10^{-8}$ eV cm, ½; cubic 123.7 meV, $4.65\times10^{-8}$, ½ | 2 | Kagamihara et al. 2004 via review Table XXV |
| Al | 220 meV, $2.32\times10^{-8}$ eV cm, 1 | 4 | Achatz et al. 2008 via review Table XXIV |

$N_D^+(n)$ and $N_A^-(p)$ replace the static net doping in every Poisson (devsim node models
`IonizedDonors` / `IonizedAcceptors`, the batched and 2D engines) and in the small-signal
capacitance (3.5). The contacts use the neutral state (`neutral_state`: bisection on
$p-n+N_D^+-N_A^- = 0$ with $np=n_i^2$). Assumed: $g$ independent of $T$, no compensation
term in $N$, no impurity band, quasi-static occupancy. Not carried by the dark-current,
FLASH and fluence engines (those run on devsim with it) and exclusive with explicit radiation
traps. **Status:** parameters **literature values**; implementation **verified** across
engines (`tests/test_incomplete_ionization_engines.py`: $C$ to $3.2\times10^{-8}$, CCE to
$10^{-9}$); the physical effect **predicted**. It is not the default because both measured-C-V
presets were fitted without it. Effect on the default device (USAGE §5.3): $V_{bi}$ 3.198 to
3.142 V, $\varepsilon A/C(0)$ 1.655 to 1.619 µm.

**Worked example** (`worked_examples.py device`, sections 2-3): legacy $f_A = 0.132222$, by hand
$(1-0.8808)\times8.554\times10^{-4} + 0.8808\times0.150$; opt-in, neutral, 300 K:
$\Delta E_{\rm Al}(10^{19}) = 0.17002$ eV and $N_A^-/N_A = 0.029131$ (code and the hand
quadratic $p^2+p_1p-p_1N_A = 0$ agree), $N_D^+/N_D$ = 0.95875 / 0.86582 / 0.51178 at
$8.19\times10^{15}$ / $4.52\times10^{16}$ / $10^{19}$ (code and a hand two-site bisection agree
to every digit).

### 3.5 Capacitance, depletion width and the C-V inversion

- **Definition used by the facades.** $W_C = \varepsilon A/C$ with the *small-signal*
  $C = dQ/dV$, $Q$ the terminal hole charge on the anode side
  (`etna/core/cv_analysis.py` `small_signal_capacitance`):
  $C = q\sum_iV_i\,[\,p_i(V_c-\delta)-p_i(V_c+\delta)\,]/(2\delta)$, $\delta = 1$ mV, two
  DC solves (with incomplete ionization the summed quantity is $p-N_A^-$). The engine
  computes the same derivative exactly from one linear solve with the Newton Jacobian
  (`PoissonBatch.small_signal_capacitance`): $C = (qN_0/V_t)\sum_iv_i\hat p_i\,du_i/du_c$.
  This is the low-frequency (quasi-static) C-V.
  **Docstring note:** the `cv_sweep` and module docstrings of `cv_analysis.py` say "devsim
  AC solve"; the code uses the ±1 mV central difference (AC is avoided because devsim's
  `delete_circuit` corrupts state, as the function's own docstring explains).
- **Depletion edge** (`extract_depletion_width_numerical`): the first n-side node with
  $n/N_D\ge0.5$, interpolated; $W_C\ge W_{\rm edge}$ (9.5185 vs 9.4803 µm at $-10$ V).
- **Analytic forms** (`etna/core/analytical.py`): one-sided abrupt
  $W = \sqrt{2\varepsilon(V_{bi}-V)/(qN_D)}$ (`depletion_width`, clamped to $t_{\rm epi}$),
  and the graded moment integral
  $V_{fd} = \frac{q}{\varepsilon}\int_0^{t}y\,N_D(y)\,dy - V_{bi}$
  (`full_depletion_voltage_graded`, closed form
  $\int_0^ty N\,dy = N_bt^2/2 + (N_j-N_b)\,[L^2-e^{-t/L}(Lt+L^2)]$; trapezoid version
  `full_depletion_voltage_profile`). **Docstring note:** `depletion_width` quotes
  "W(0 V) ~1.7 µm", true only for the old uniform calibration ($10^{15}$); on the default
  graded stack the uniform formula gives 0.63 µm ($N_j$) or 9.77 µm ($N_b$). The helper
  `poisson.extract_depletion_width` (analytic $W$ on the obsolete scalar `N_D`) has no caller
  and must not be used on the default stack.
- **Doping from C-V** (`etna/measurements/cv.py`). (a) Differential profile
  $N(W) = 2/\big(q\varepsilon\,d(1/C'^2)/dV_{\rm rev}\big)$, $C' = C/A$, $W = \varepsilon/C'$,
  local polynomial slope (`depletion_profile`); it smears over about one Debye length
  $L_D = \sqrt{\varepsilon k_BT/(q^2N)}$ and diverges at full depletion. (b) Model inversion
  (`invert_cv`): least squares on $\log_{10}N_D$ at knots with a second-difference
  regulariser, every trial a full batched Poisson plus small-signal C. The tabulated preset
  is that inversion; the parametric default was fitted the same way (forward model), not by
  Mott-Schottky.

**Status:** Poisson, contacts and $C$ **verified** code to code (engine vs devsim
$|\Delta C|/C\lesssim3\times10^{-8}$, `tests/test_cv_small_signal.py`; parity of the transition
and table profiles $|\Delta W|/W\le1.9\times10^{-15}$); the solver scale **cross-checked** on a
Schottky diode of another wafer (Bruzzi & Verroi 2023 [BV23], $W$ +3.1 % at $-40$ V, +1.0 %
with $V_{bi}$ matched; `scripts/validation/cv_vs_bruzzi2023.py`); the C-V of this device
**calibrated**.

**Worked example** (`worked_examples.py device`, sections 5-6). The table compares the full
Poisson solution with the analytic forms:

| $V$ (V) | $\varepsilon A/C$ (µm) | $W_{\rm edge}$ | uniform $N_j$ | uniform $N_b$ | graded moment |
| --- | --- | --- | --- | --- | --- |
| 0 | 1.6547 | 1.6493 | 0.6283 | 9.7746 | 1.7564 |
| -0.5 | 2.5402 | 2.4691 | 0.6784 | 10.000 | 2.9855 |
| -2 | 6.9390 | 6.7922 | 0.8104 | 10.000 | 7.4314 |
| -10 | 9.5185 | 9.4803 | 1.3056 | 10.000 | 9.5525 |

The uniform formula is useless on the graded profile; the moment inversion is within 0.1 µm at
0 and $-10$ V but overshoots in the knee, where Debye tails ($L_D$ 0.4-0.6 µm) matter. The
differential profile of the computed C-V recovers the bulk $3.38\times10^{13}$ to 1-7 % between
$-1$ and $-2$ V ($3.80\times10^{13}$ at $W$ = 4.17 µm against the true $3.48\times10^{13}$;
$3.62\times10^{13}$ at 6.94 µm against $3.383\times10^{13}$) and fails by $\times2$ within one
$L_D$ of the transition (8.5 µm), the expected Debye smearing.

## 4. Drift-diffusion, Scharfetter-Gummel and the box method

### 4.1 Continuous equations

Together with the Poisson equation of chapter 3, ETNA solves the electron and hole continuity
equations

$$
\frac{\partial n}{\partial t} = \frac{1}{q}\nabla\cdot\mathbf J_n - U + G,\qquad
\frac{\partial p}{\partial t} = -\frac{1}{q}\nabla\cdot\mathbf J_p - U + G,
$$

$$
\mathbf J_n = q\mu_n n\,\mathbf E + qD_n\nabla n,\qquad
\mathbf J_p = q\mu_p p\,\mathbf E - qD_p\nabla p,\qquad D = \mu V_t ,
$$

with the Einstein relation of non-degenerate statistics, $U$ the net recombination of
chapter 6 and $G$ the radiation (or effective) generation. There is no energy balance or
hydrodynamic transport, no thermodiffusion, and no band-gap narrowing in the solved
equations.

| solver | files | unknowns, coupling |
| --- | --- | --- |
| devsim (reference) | `etna/core/drift_diffusion.py` `setup_sic_drift_diffusion` | $\psi,n,p$, fully coupled Newton |
| batched 1D (numba) | `etna/core/batched/poisson.py`, `transport.py`, `models_dd.py` | Gummel: Poisson with nodal quasi-Fermi levels, then $n$, then $p$; with the opt-in models a full Newton on $(n,p)$ at frozen $\psi$ |
| 2D laterally uniform (numba) | `etna/core/engine2d/` | dark state by quasi-Fermi Poisson, illuminated state by Gummel with the dark LUs as chord corrections |
| general 2D (numba, opt-in) | `etna/core/fv2d/` | as engine2d on any tensor mesh, Cartesian or $(r,z)$, any contacts |
| track engine (numba, opt-in) | `etna/core/track/` | $(\delta\psi,n,p)$ or the excess $(\delta\psi,\delta n,\delta p)$, fully coupled Newton, TR-BDF2 |

All of them use the discretisation below, which is why the engines reproduce devsim to
round-off inside their checked boxes (16.6).

### 4.2 The box (finite-volume) method

Each node $k$ owns a control volume $V_k$ (the dual cell). In 1D, with $h_i = x_{i+1}-x_i$,
$V_i = (h_{i-1}+h_i)/2$, half cells at the contacts. On ETNA's 2D tensor meshes devsim's dual
cell is the rectangle $\Delta x^{\rm dual}_i\Delta y^{\rm dual}_j$; this is asserted, not
assumed, in `etna/core/engine2d/grid.py` (devsim's diagonal edges carry a relative coupling
of about $10^{-33}$). Each edge $k\to l$ carries the geometric factor $S_{kl}/h_{kl}$
(devsim `EdgeCouple`/`EdgeLength`). In $(r,z)$ (`etna/core/track/grid.py`,
`etna/core/fv2d/mesh.py` `box_geometry`) the volumes are
$V_{ij} = \pi(r_{f,i+1/2}^2-r_{f,i-1/2}^2)\Delta z_j$, the radial faces $2\pi r_{f,i+1/2}\Delta z_j$.

Integrating the electron equation over $V_k$:

$$
V_k\,\frac{dn_k}{dt} = \sum_{l\in\mathcal N(k)}\frac{S_{kl}}{q}\,J_{n,kl} + V_k\,(G_k-U_k),
$$

so sources are lumped at nodes ("node models") and fluxes live on edges ("edge models").
Two consequences recur in this manual:

- the quadrature consistent with the solver is the **box rule** $\int G\,dV\approx\sum_kV_kG_k$
  (`integrate_nodes_box_2d` in `etna/core/charge_collection_2d.py`; in 1D it equals the
  trapezoid rule that `compute_cce_from_dd` uses);
- a contact node has a Dirichlet row (3.2) that *replaces* its continuity equation, so
  generation assigned to it is never injected (`injected_generation`), and generation in the
  contact half-cell is lumped onto the contact. This is the origin of the contact-cell
  finding of 16.5.

### 4.3 Scharfetter-Gummel edge currents

On an edge from node $k$ (`n0`) to node $l$ (`n1`) of length $h$, with
$\Delta = (\psi_k-\psi_l)/V_t$ and $B(x) = x/(e^x-1)$, devsim's
`simple_dd.CreateElectronCurrent` / `CreateHoleCurrent` (used unchanged by
`setup_sic_drift_diffusion`) are

$$
J_{n,kl} = \frac{q\mu_nV_t}{h}\big[n_l\,B(-\Delta) - n_k\,B(\Delta)\big],\qquad
J_{p,kl} = \frac{q\mu_pV_t}{h}\big[p_k\,B(-\Delta) - p_l\,B(\Delta)\big],
$$

positive along $k\to l$. devsim evaluates the electron bracket as
`kahan3(n1*Bern01, n1*vdiff, -n0*Bern01)`, i.e. with $B(-x) = B(x)+x$ and compensated
summation. The scheme is exact for constant current, field and mobility along the edge; it
tends to central differences as $\Delta\to0$ and to upwinded drift for $|\Delta|\gg1$. The
edge mobility is the mean of the node values, $\mu_{kl} = (\mu_k+\mu_l)/2$
(`add_local_mobility_models` in `etna/core/device.py`; the same in the numba kernels).

The numba Bernoulli function (`bernoulli` / `_bern`, `etna/core/batched/transport.py`) is
$1-x/2$ for $|x|<10^{-8}$, $x/\mathrm{expm1}(x)$ for $|x|\le700$, $0$ above and $-x$ below.
The track engine's excess unknowns also need $B(x+d)-B(x)$ without cancellation
(`bern_diff`, series and `expm1` forms, $10^{-14}$ relative against a 50-digit reference,
USAGE §19.1).

### 4.4 Discrete equations as assembled

The batched 1D rows (`etna/core/batched/transport.py`), $v_i$ the node volume:

$$
J_{n,i}-J_{n,i-1} = q\,v_i\,(U_i-G_i),\qquad J_{p,i}-J_{p,i-1} = -q\,v_i\,(U_i-G_i).
$$

The 2D kernels (`etna/core/engine2d/kernels.py`, `etna/core/fv2d/kernels.py`) divide every
equation by $q$ and use $u=\psi/V_t$:

$$
F^\psi_k = \lambda\sum_bg_{kb}(u_k-u_b) - V_k(p-n+D)_k,\qquad \lambda = \frac{\varepsilon V_t}{q},\ g = \frac{S}{h},
$$

$$
R^n_k = \sum_bk_n\big[n_bB(u_b-u_k) - n_kB(u_k-u_b)\big] - V_kU_k + V_kG_k,\qquad k_n = \mu_nV_tg,
$$

holes symmetric, with $n = e^{u-u_n+\ln n_i}$, $p = e^{u_p-u+\ln n_i}$ in the quasi-Fermi
form. Unless a model replaces it, $U$ is **midgap SRH**, $n_1 = p_1 = n_i$ (devsim `USRH`, all
engines). The dark-current model (6.3) replaces it by a Z1/2-level SRH plus an effective
generation; the FLASH and track paths add Auger (6.2).

### 4.5 Validity

Local drift-diffusion with Boltzmann statistics: no hot-carrier, ballistic or degenerate
effects, and no tunnelling other than the Hurkx factor of 6.3. In the sub-nm cores of
heavy-ion tracks the first-cell density can exceed the atomic density (USAGE §19); the track
results state that drift-diffusion is outside its validity there. **Status:** the
discretisation is devsim's; the engines are **verified** against it (16.6).

**Worked example** (`worked_examples.py sg`). The code's $B(x)$ equals $x/(e^x-1)$ at every
printed digit from $x=-50$ to $50$ ($B(1) = 0.5819767068693263$,
$B(-1) = 1.581976706869326$, $B(50) = 9.643749239819590\times10^{-21}$) and
$B(-x)-B(x)-x\le6\times10^{-17}$. For uniform $n = 10^{14}$ cm$^{-3}$, 0.01 V across 0.1 µm
and $\mu_n = 1000$ cm$^2$/Vs the SG flux is $16.02177$ A/cm$^2$, equal to the drift current
$q\mu_nnE$: no numerical diffusion for this exact solution.

## 5. Mobility

The parameter sets live in registries in `etna/core/physics_models.py`; each carries a
`Provenance` record stating what was read and where ("verified-primary" = numbers read in
the primary source, "verified-secondary" = read in a review that restates it, "assumption").
`PhysicsModels()` is the default physics: the device's own low-field set, no velocity
saturation, no anisotropy model, no avalanche, no band-gap narrowing.

### 5.1 Low-field, doping- and temperature-dependent mobility

Two forms, with $N$ the **total** local doping $N_D + N_A$ and $t = T/300$
(`LowFieldMobility.mobility`):

$$
\text{Caughey-Thomas ("ct")}:\quad \mu = \mu_{\min}t^{\beta} + \frac{\mu_{\max}t^{\gamma}-\mu_{\min}t^{\beta}}{1+(N/N_{\rm ref})^{\alpha}},
$$

$$
\text{Arora-type ("arora")}:\quad \mu = \frac{(\mu_{\max}-\mu_{\min})\,t^{\gamma}}{1+(N/N_{\rm ref})^{\alpha}\,t^{\delta}} + \mu_{\min}t^{\beta}.
$$

The solvers evaluate it per node at the local total doping (the total $N_A$ is the ionized
node value divided by $f_A$) and average it on edges (4.3)
(`add_local_mobility_models`, `etna/core/device.py`; `mobility_caughey_thomas_T`,
`etna/core/sic_material.py`). `DeviceConfig(mobility=...)` selects the set:

| set | $\mu_{\max},\mu_{\min},N_{\rm ref},\alpha,\gamma,\beta$ electrons | holes | direction | provenance |
| --- | --- | --- | --- | --- |
| `c_axis` (**default** since 2026-10-21) | 1190, 10, $10^{17}$, 0.54, -2.40, 0 | 125, 15.9, $1.76\times10^{19}$, 0.34, -2.15, 0 | parallel to c | electrons: Ishikawa 2025 thesis Eq. (2.16), drift mobility parallel to c [ISH25]; $\gamma_n$ is Schaffer's $\langle0001\rangle$ $-2.40$, not the thesis's $-3.4$; holes: the legacy set |
| `legacy_perp` | 950, 40, $1.94\times10^{17}$, 0.61, -2.40, 0 | as above | perpendicular to c | Schaffer et al. 1994 via Ayalew Table 3.5 [AYA] |
| `c_axis_holes_0p9` (opt-in) | as `c_axis` | holes $\times0.90$ (112.5, 14.31, ...) | parallel | thesis drift ratio [ISH25] |
| `ishikawa2025` (opt-in) | 1190, 10, $10^{17}$, 0.54, -3.4, 0.8; $\delta$ -2.2 | 83, 20, $6.4\times10^{18}$, 0.83, -2.3, -2.2; $\delta$ -0.9 | parallel | full thesis laws, fitted 300-600 K (warns outside) |
| `ayalew2004` (opt-in) | 950, 40, $1.94\times10^{17}$, 0.61, -2.40, -0.5 | 125, 15.9, $1.76\times10^{19}$, 0.34, -2.15, -0.5 | perpendicular | Ayalew Table 3.5, eqs. 3.77-3.79 |

**Why parallel to c.** The 1D diode is a vertical device on a (0001) wafer, so its field and
current are along c (the 4° off-cut is neglected). Until 2026-10-21 ETNA used the
perpendicular (Hall, basal-plane) set for vertical transport; the `c_axis` electrons are
+24 % at the epi bulk. The provenance notes record an uncertainty of about ±8 % for the
electrons (Hall factor) and $-35/+0$ % along c for the holes, whose magnitude and anisotropy
disagree between the primaries; electrons in the p$^+$ (Al-doped) are an extrapolation of an
N-doped fit.

**Deviation.** In `legacy_perp`, $\mu_{\min}$ is held T-independent (Ayalew: $(T/300)^{-0.5}$);
the hole $\gamma = -2.15$ is borrowed from the electrons.

**2D.** The 2D solvers use the same isotropic node mobility on lateral (basal-plane) edges:
lateral anisotropy is not modelled.

### 5.2 Anisotropy (opt-in)

`ANISOTROPY` gives $\mu_\parallel/\mu_\perp$ per carrier: `ayalew2004` (1.2, 1.0; Ayalew eqs.
3.80-3.81) and `burin2024` ($1/0.83$, $1/1.15$; review Table XXVII [BUR24r]). It multiplies a
*perpendicular* low-field set to obtain vertical (parallel) transport; `PhysicsModels` refuses
it on a set that is already parallel to c.

### 5.3 High-field mobility (opt-in)

Caughey-Thomas / Canali form on every edge, evaluated at the edge field
$E_e = |\psi_i-\psi_{i+1}|/h_e$ (`HighFieldMobility.mobility`; numba
`edge_mobility_kernel` in `etna/core/batched/models_dd.py`; devsim twin in
`etna/core/devsim_models.py`, where the field factor is inlined into the SG current so that
its derivatives enter the Jacobian):

$$
\mu(E) = \frac{\mu_0}{\big[1+(\mu_0E/v_{\rm sat})^{\beta}\big]^{1/\beta}},\qquad
v_{\rm sat}(T) = v_{\rm sat}\,t^{\delta},\quad \beta(T) = \beta\,t^{\gamma_\beta}.
$$

| set | electrons $v_{\rm sat}$ perp / par (cm/s), $\beta$, $\delta$, $\gamma_\beta$ | holes | provenance |
| --- | --- | --- | --- |
| `lades2000` | $2.2\times10^7$ / $1.64\times10^7$, 1.2, -0.44, 1.0 | $1.1\times10^7$ / $6.5\times10^6$, 1.2 (**assumption**: electron $\beta$, exponents) | Lades 2000 via Ayalew Table 3.6; parallel and hole $v_{\rm sat}$ Nilsson 2000 MC via review Tables XXXIII-XXXIV |
| `ishikawa2025` | $1.7\times10^7$ both, 0.76, no T law | the `lades2000` holes | thesis Table 3.1 (pulsed conductance) |

The parallel $v_{\rm sat}$ is used when the vertical transport is parallel to c (the default
set, or a perpendicular set with an anisotropy model). Velocity saturation is used by the
opt-in facades `run_avalanche`, `run_transit_time` (`etna/api/models.py`) and as a
frozen-field bracket in the track cross-sections (`mu_field`); **no default facade has it**.
At the operating fields of the default device ($4\times10^4$ V/cm average at $-40$ V,
$1.7\times10^5$ V/cm at the junction peak) the electron $\mu_0E/v_{\rm sat}$ is 3 to 12, so the
drift velocities of the default solves are overestimated. Collected charge at low injection
should be insensitive to this (the carriers still reach the contacts; not computed here); transit times, pulse shapes
and the plasma erosion of tracks are not.

### 5.4 Band-gap narrowing (tabulated, not coupled)

`BANDGAP_NARROWING["lades2000"]`: Slotboom $\Delta E_g = C(F+\sqrt{F^2+G})$,
$F = \ln(N/10^{17})$, $C_n = 2.0\times10^{-2}$, $C_p = 9.0\times10^{-3}$ eV, $G = 0.5$ (Lades 2000,
fit to Lindefelt 1998, via Ayalew eq. 3.72 / Table 3.4). Selecting it in a solver raises
`NotImplementedError` (`etna/core/devsim_models.py`).

**Status:** parameter sets **literature values**; engine and devsim twin **verified**
against each other with the opt-in models (`tests/test_physics_models.py`); the choice of the
c-axis set is a documented physics decision whose effect on CCE is **predicted**
([USAGE §18](USAGE.md#18-physical-model-library-opt-in-high-field-mobility-anisotropy-avalanche-and-breakdown)).

**Worked example** (`worked_examples.py materials`): at 300 K the default electron mobility is
1174.44 / 600.00 / 100.61 cm$^2$/Vs at $N$ = $3.38\times10^{13}$ / $10^{17}$ / $10^{19}$, equal to
$10 + 1180/(1+(N/10^{17})^{0.54})$ by hand; `legacy_perp` gives 945.39 / 585.73 / 115.35, the
holes 123.77 / 108.96 / 75.68. At 350 K the bulk electron mobility drops to 811.3 cm$^2$/Vs
($\mu_{\max}\,(350/300)^{-2.40}$).

## 6. Recombination, generation and the dark current

### 6.1 Shockley-Read-Hall

$$
U_{\rm SRH} = \frac{np-n_i^2}{\tau_p(n+n_1)+\tau_n(p+p_1)},
$$

midgap $n_1 = p_1 = n_i$ in every transport solve (devsim `USRH`,
`setup_sic_drift_diffusion`; the numba kernels), with the scalar lifetimes
$\tau_n = 1$ ns, $\tau_p = 600$ ns at 300 K (Ioffe) and $\tau(T) = \tau_{300}(T/300)^{1.72}$
(`srh_lifetime`). The lifetimes are region-wide scalars, the same in the p$^+$, epi and n$^+$.
No doping dependence (Scharfetter) is applied. Radiation damage shortens them (12.6).

With $n_i\approx5\times10^{-9}$ cm$^{-3}$, midgap SRH generation in the depletion region is
$\sim qn_iW/(2\tau)\approx10^{-49}$ A for this device (module docstring of
`etna/core/dark_current.py`): it predicts no measurable leakage, which is why the dark current
uses the effective model of 6.3.

**Status:** **literature values** (Ioffe); the lifetimes of this epi are not measured. The
tilted-alpha study of USAGE §19.3 needs about 23 ns (band 17-42 ns) for $\tau_n+\tau_p$ to fit
the shape of one microdosimeter curve, a **calibration** that the defaults do not adopt.

### 6.2 Auger

$$
U_{\rm Aug} = (C_nn + C_pp)\,(np-n_i^2),\qquad C_n = 5\times10^{-31},\ C_p = 2\times10^{-31}\ \text{cm}^6/\text{s},
$$

(`add_auger_recombination`, `etna/core/flash_recombination.py`; numba
`etna/core/batched/dark.py`; track `kernels.py`). Galeckas et al. 1997 measured them by
free-carrier absorption at 300 K, $\Delta N = 2\times10^{18}$-$2\times10^{19}$ cm$^{-3}$
[GAL97]; later measurements span $10^{-31}$-$2\times10^{-30}$ with injection and temperature
dependences (review Table XXII [BUR24r]). ETNA uses 300 K constants at every $T$. Auger is
**off** in the default CCE and C-V solves; it is on in the FLASH study (chapter 11) and the
track engine, which also offers the injection-dependent reduction
$f = \min(1, A\Delta^{-B}/(C_n+C_p))$, $\Delta = \min(n,p)$ (`auger_model="injection"`, Tanaka
2023 via the review). Radiative recombination is not modelled.

**Worked example** (`worked_examples.py recombination`): in the $10^{19}$ n$^+$ the minority-hole
Auger lifetime $1/(C_nN_D^2) = 2.0\times10^{-8}$ s is shorter than $\tau_p = 600$ ns: this is the
dose-rate-independent CCE offset of chapter 11.

### 6.3 The dark-current model ("TAT", calibrated)

`setup_tat_model` (`etna/core/dark_current.py`; numba twin `etna/core/batched/dark.py`)
replaces $U$ by

$$
U_{\rm TAT} = \frac{np-n_i^2}{\dfrac{\tau_p}{1+\Gamma}(n+n_1^{t}) + \dfrac{\tau_n}{1+\Gamma}(p+p_1^{t})}
\;-\; G_0\,(1+\Gamma)\,\min\!\left(\frac{E_{\rm node}}{E_{\rm ref}},1\right)\,\mathrm{step}(x_{\rm end}-x),
$$

with

- the Z1/2 trap level $E_t = E_c-0.67$ eV: $n_1^t = n_i\,e^{(E_g/2-E_t)/k_BT}$,
  $p_1^t = n_i\,e^{-(E_g/2-E_t)/k_BT}$ (the intrinsic level taken at midgap);
- the Hurkx et al. 1992 field enhancement (cited in the code; the paper was not opened, the
  closed form is derived in the `hurkx_gamma` docstring), low-field closed form

$$
\Gamma = 2\sqrt{3\pi}\,\frac{|E|}{F_\Gamma}\exp\!\left[\left(\frac{E}{F_\Gamma}\right)^2\right],\qquad
F_\Gamma = \frac{\sqrt{24\,m_t\,(k_BT)^3}}{q\hbar},
$$

  $m_t = 0.25\,m_0$ (**assumption**), exponent capped at 700 (`hurkx_gamma`), valid while the
  integrand's saddle point $u^* = (2E_t/3k_BTK)^2\le1$, i.e. $E\le E_{\rm hurkx,max}$
  ($1.08\times10^6$ V/cm at 300 K); $\Gamma_n = \Gamma_p$ (one tunnelling mass);
- the effective generation $G_0 = N_t\,n_i(T)/n_i(300)$, $N_t = 1.355575\times10^{12}$
  cm$^{-3}$ s$^{-1}$, $E_{\rm ref} = 3.5\times10^4$ V/cm, and $E_{\rm node}$ the node field
  (end nodes: their edge; interior: the mean of the two edges; floored at $10^3$ V/cm). The
  devsim path iterates the frozen node field to self-consistency; the engine re-evaluates it at
  every Gummel iteration, which converges to the same fixed point. The step confines the
  effective generation to the p$^+$ and epi when an n$^+$ layer exists.

**What is fitted and what is assumed.** $N_t$ is **calibrated** so that
$I(300\text{ K}, -60\text{ V}, 4\text{ mm}^2) = 18.000$ pA, the measured *upper bound*
"< 18 pA up to $-60$ V" used as an equality. The bias shape comes from the heuristic
$f_{\rm depl} = \min(E/E_{\rm ref},1)$ (**assumption**: $E_{\rm ref}$ is half the ~70 kV/cm
built-in junction field) times $(1+\Gamma(E))$. The temperature law $N_t\propto n_i(T)$
(activation $E_g/2$) is an **assumption**; with $\Gamma(T)$ the model gives $\times15.1$ from 300
to 313 K at $-30$ V. The SRH/TAT/surface split is bookkeeping, not physics. Leakage scales with
the area passed in.

**Deviation.** The closed form for $\Gamma$ is the saddle-point evaluation of Hurkx's integral
$\Gamma = a\int_0^1e^{au-Ku^{3/2}}du$; at these fields it is 1.13-1.32$\times$ the integral with
$\Delta E = E_t$ (`tests/test_dark_current.py`). The absolute scale goes into the $N_t$ fit; the
bias shape keeps that $\approx30$ % approximation.

**Worked example** (`worked_examples.py recombination`): $F_\Gamma = 3.6878\times10^5$ V/cm in
the code, $3.6884\times10^5$ by hand with CODATA constants (the code's rounded $\hbar$, $q$, $m_0$
make the 1.6e-4 difference); $\Gamma$ = 0.848 / 1.792 / 3.501 at $5\times10^4$ / $10^5$ /
$1.7\times10^5$ V/cm (code = hand); $n_1^t = 6.71\times10^7$, $p_1^t = 3.7\times10^{-25}$ cm$^{-3}$.

### 6.4 Hecht (analytic CCE reference)

`hecht_cce` (`etna/core/charge_collection.py`): two-carrier Hecht for uniform generation in
a uniform field $E = |V|/d$, averaged over the generation depth so that each pair counts once
(audit C1):

$$
\mathrm{CCE} = \sum_{c=e,h}\left[\frac{\lambda_c}{d} - \left(\frac{\lambda_c}{d}\right)^2\left(1-e^{-d/\lambda_c}\right)\right],\qquad \lambda_c = \mu_c\tau_c\,\frac{|V|}{d}.
$$

Each term tends to 1/2 for $\lambda\gg d$. `hecht_cce_partial_depletion` adds a diffusion term
$(L/t)(1-e^{-t/L})$ for the undepleted thickness $t$ with a hole diffusion length
$L = 7$ µm (a literature CCE fit recorded only by its NIM A PII, authors not recorded).

**Docstring note.** The `hecht_cce` docstring quotes default mobilities 945.4 / 123.7
cm$^2$/Vs (the `legacy_perp` set); the code calls `mobility_caughey_thomas_T` with the
default `c_axis` set and uses 1174.44 / 123.77.

**Worked example** (`worked_examples.py recombination`): $V = 40$ V, $d = 10$ µm:
$\lambda_e = 470$ µm, $\lambda_h = 29705$ µm, CCE = 0.996415 in the code and by hand.

## 7. Impact ionisation and breakdown (opt-in)

### 7.1 Coefficients

All sets have the form (`AlphaParams`, `etna/core/physics_models.py`; numba `_alpha`,
`etna/core/avalanche.py`)

$$
\alpha(E,T) = a\,[1+c(T-300)]\,\exp\!\left\{-\left[\frac{b\,[1+d(T-300)]}{|E|}\right]^m\right\},
$$

set to 0 below $|E| = 1$ V/cm, with $E$ the edge field. The 1D solvers use the
parallel-to-c coefficients and refuse a set without them; a set used outside its measured
temperature range needs `allow_T_extrapolation=True`.

| set | electrons $a$ (1/cm), $b$ (V/cm), $m$ | holes | range | provenance |
| --- | --- | --- | --- | --- |
| `hatakeyama2004` (default of `run_breakdown`) | $1.76\times10^8$, $3.30\times10^7$, 1 | $3.41\times10^8$, $2.50\times10^7$, 1 | 300 K, 2-5 MV/cm | Hatakeyama et al., APL 85, 1380 (2004) Table I, primary [HAT04]; also perpendicular: $2.10\times10^7$/$1.70\times10^7$, $2.96\times10^7$/$1.60\times10^7$ |
| `niwa2015` | $0.143\times10^6$, $4.93\times10^6$, 2.37 | $3.14\times10^6$, $11.8\times10^6$, 1.02, $c$ 6.3e-3, $d$ 1.23e-3 | 300-448 K (assumed) | Niwa et al. 2015 via review Table XIII |
| `zhao2019` | $0.339\times10^6$, $5.15\times10^6$, 2.37 | $3.56\times10^6$, $11.7\times10^6$, 1.02, 6.19e-3, 1.15e-3 | 150-673 K | via review |
| `loh2008` | $2.78\times10^6$, $10.5\times10^6$, 1.37 | $3.51\times10^6$, $10.3\times10^6$, 1.09 | 300 K | via review (direction assumed) |
| `stefanakis2021` | $2.8\times10^6$, $20.7\times10^6$, 1, -1.0e-3, -0.29e-3 | $2.5\times10^6$, $12.1\times10^6$, 1, 1.74e-3, 0.59e-3 | 300-573 K (assumed) | via review |
| `ayalew2004_legacy` | $3.44\times10^6$, $2.58\times10^7$ | $3.5\times10^6$, $1.7\times10^7$ | 300 K | **legacy, not recommended** (6H-derived $b_n$) |

Hatakeyama's anisotropic interpolation for arbitrary field direction
(`hatakeyama_anisotropic_alpha`) is implemented for 2D/3D use but not coupled into a solver.

### 7.2 McIntyre ionisation integrals and breakdown

With $\alpha_{n,p}$ constant on each mesh edge and
$\varphi(x) = \int_0^x(\alpha_n-\alpha_p)\,dx'$ (anode at 0, electrons drifting to $+x$),
`etna/core/avalanche.py` evaluates the McIntyre (1966) relations

$$
I_n = \int_0^W\alpha_ne^{-\varphi}\,dx,\qquad I_p = \int_0^W\alpha_pe^{\varphi(W)-\varphi}\,dx,
$$

$$
M(x) = \frac{e^{\varphi(W)-\varphi(x)}}{1-I_p},\qquad M_n = M(0) = \frac{1}{1-I_n},\quad M_p = M(W) = \frac{1}{1-I_p},
$$

with each edge integrated in closed form: on an edge of length $h$ with
$z = (\alpha_n-\alpha_p)h$ the factor is $h\,(1-e^{-z})/z$ (`_f`, series for $|z|<10^{-6}$).
The identity $1-I_n = e^{-\varphi(W)}(1-I_p)$ holds to round-off, so both integrals reach 1 at
the same bias, the breakdown voltage $V_{BD}$. A distributed source has the
generation-weighted gain $M_G = \int GM\,dx/\int G\,dx$.

`breakdown_voltage` marches each device in reverse bias with batched Poisson solves (sub-steps
of at most 10 V, scaled down below 300 K), predicts the next bias by a secant on $\ln I_p$, and
closes the bracket by Illinois regula falsi to $10^{-3}$ V. Two error estimates are reported:
`quadrature_error_V` ($|\ln I(h)-\ln I(2h)|/3$ over the local slope; quadrature only, about
100$\times$ below the real mesh error) and, in `run_breakdown`, `mesh_error_V`
$= |V_{BD}(h)-V_{BD}(h/2)|$ from a full re-solve on the bisected mesh (convergence about first
order, so the base value is within about $2\times$ that number of the converged one).

In the batched DD solve with avalanche (`etna/core/batched/models_dd.py`) the impact
generation $g_e = \alpha_n|J_{n,e}| + \alpha_p|J_{p,e}|$ (A/cm$^3$) is integrated over the edge
volume (node $i$ receives $\tfrac12(h_{i-1}g_{i-1}+h_ig_i)$) and solved by a full Newton on
$(n,p)$ at frozen $\psi$ (the avalanche feedback is not lagged; a Gummel-lagged term converges
with ratio $1-1/M$).

**Validity and status.** Parallel-plane 1D: no edge termination, so $V_{BD}$ of a real device
is lower; the default device's value is an **upper bound** and a **prediction**. The
coefficients are **literature values**; against Hatakeyama's measured
$V_{BD} = 1940\,(10^{16}/N_D)^{0.8}$ V on the mesas the coefficients were fitted on, the model
is consistent (not validated); the other sets show the model spread
([USAGE §18](USAGE.md#18-physical-model-library-opt-in-high-field-mobility-anisotropy-avalanche-and-breakdown)).
The engine vs devsim twin is **verified** (`tests/test_physics_models.py`,
`tests/test_physics_models_devsim.py`).

**Worked example** (`worked_examples.py avalanche`): `hatakeyama2004` along c gives
$\alpha_n = 12.01$ cm$^{-1}$ and $\alpha_p = 1271$ cm$^{-1}$ at 2 MV/cm, equal to
$1.76\times10^8e^{-3.30\times10^7/E}$ and $3.41\times10^8e^{-2.50\times10^7/E}$; the default
dosimeter breaks down (parallel plane) at $V_{BD} = 2225.78$ V with `mesh_error_V` 0.105 V and
a peak field of $2.29\times10^6$ V/cm.

## 8. Generation by radiation

### 8.1 Pair creation energy and dose rate

Every deposit is converted to pairs with $E_{\rm pair} = 7.8$ eV (`SiC4H_Parameters.E_pair_eV`,
"measured W-value (Lebedev; ~7.7-7.8 eV)"; no full reference is recorded in the code, so the
value is a **literature value, source not verified here**). The Allpix Squared default for
4H-SiC, 7.6 eV, is quoted in `etna/core/mc_response.py` for comparison. A uniform dose rate
$\dot D$ in SiC gives (`dose_rate_to_generation`, `etna/core/generation_profiles.py`)

$$
G = \frac{\dot D\,\rho\times10^4}{E_{\rm pair}\times1.602\times10^{-12}}\quad[\text{cm}^{-3}\,\text{s}^{-1}],
$$

with $\dot D$ in Gy/s, $\rho = 3.21$ g/cm$^3$, $10^4$ erg/g per Gy and $1.602\times10^{-12}$ erg/eV.

### 8.2 NIST stopping powers and SiC by Bragg additivity

`etna/core/stopping_power.py` reads NIST PSTAR (protons) and ASTAR (alphas) tables for Si,
amorphous carbon, graphite, liquid water and ICRU muscle, downloaded 2026-10-20 with URLs,
request parameters and SHA-256 recorded in `data/nist_star/README.md` [NIST]. NIST has no SiC
table; SiC is composed by Bragg additivity of the mass stopping powers,

$$
\left(\frac{S}{\rho}\right)_{\rm SiC} = w_{\rm Si}\left(\frac{S}{\rho}\right)_{\rm Si} + w_{\rm C}\left(\frac{S}{\rho}\right)_{\rm C},\qquad w_{\rm Si} = 0.7004,\ w_{\rm C} = 0.2996,
$$

(`sic_mass_stopping`, carbon = amorphous by default), interpolated log-log, with the linear
stopping power (LET) $L = (S/\rho)\,\rho$ (`linear_stopping`, MeV/cm; 1 MeV/cm = 0.1 keV/µm).
The CSDA range of SiC is the integral of $1/S_{\rm tot}$ of the composed table on a fine log
grid, started at 1 keV with the mass-fraction harmonic mean of the elemental NIST ranges
(`_csda_table`; the same integration reproduces NIST's own Si range to $<0.1$ %, tested).
The residual energy after a depth $s$ is read from the inverted range table
(`residual_energy_MeV`). **Assumption:** Bragg additivity ignores chemical binding and phase;
the module's error estimate (a few % above 1 MeV for protons, up to ~10 % near the stopping
maximum) is an unverified literature estimate. Graphite instead of amorphous carbon changes the
composed SiC stopping by $\le1$ % above 1 MeV.

### 8.3 Alpha particles: the Bragg curve in SiC

`alpha_bragg_curve_sic`: the ionisation along the CSDA path is the *electronic* stopping
power at the residual energy, $dE/ds = L_e(E(s))$, $E(s)$ from $R(E_0)-R(E) = s\rho$, convolved
with a Gaussian of the Bohr range straggling (`bohr_range_straggling_cm`)

$$
\sigma_R^2 = \int_0^{E_0}\frac{d\Omega_B^2/dx}{L_{\rm tot}(E)^3}\,dE,\qquad
\frac{d\Omega_B^2}{dx} = 4\pi z^2e^4N_e ,
$$

($e^2 = 1.439964\times10^{-13}$ MeV cm, $N_e$ the electron density of SiC from $\rho$ and the
mass fractions; non-relativistic, free electrons, so an upper-ish estimate). The depth profile
used by the CCE sweeps (`alpha_depth_shape` with `bragg_peak_ratio="nist"`, default) is that
curve on the unit axis $u = s/R$, and `alpha_generation_profile` normalises it so that
$\int G\,dx = E_0/E_{\rm pair}$ on the grid supplied, whatever the shape.

**Deviations, stated.** (i) The shape uses $S_e$ while the normalisation is to the full $E_0$
(the nuclear share of a 5.5 MeV alpha is small but non-zero). (ii) In `cce_vs_bias` the profile
starts at the **junction** ($x_{\rm epi} = x-x_j$), at normal incidence, and its peak value is
set to the generation rate $10^{18}$ cm$^{-3}$ s$^{-1}$ (a steady-state, low-injection proxy);
the p$^+$ dead layer and the deposit beyond the epi are not generated, so this CCE is relative to
the in-epi generation (`in_mesh_generation_fraction` reports the in-epi share for the 2D scans).
(iii) The tilted track of the microdosimeter measurement (`tilted_alpha_generation`) maps depth
$x$ to path $s = (x+d_{\rm dead})/\cos\theta$ and $G(x) = (dE/ds)(s)/(\cos\theta\,E_{\rm pair})$
with an optional SiC-equivalent dead layer; the legacy box and power-law shapes are opt-outs
labelled "not data".

### 8.4 Protons and ions

`proton_generation_profile` has three input conventions (audit-v7 B-10):

- default, `dose_medium="sic"`: $\dot D$ is the dose rate in SiC; $G$ is flat and
  energy-independent (8.1); the energy is only checked to lie in 30-150 MeV, the scope of a flat
  entrance-region profile;
- `dose_medium="water"`: $\dot D_{\rm SiC} = \dot D_{\rm w}\,(S/\rho)_{\rm SiC}/(S/\rho)_{\rm w}$
  at the beam energy (Bragg-Gray, no secondary-electron transport);
- `fluence_rate_cm2_s` $\dot\Phi$: $G(x) = \dot\Phi\,L_{e,\rm SiC}(E(x))\times10^6/E_{\rm pair}$
  with $E(x)$ the residual CSDA energy, and $G=0$ beyond the range; valid at any PSTAR energy.

For a single ion the line density of pairs is
$N_L = L/E_{\rm pair}$, i.e. $L[\text{keV/µm}]\times10^7/E_{\rm pair}[\text{eV}]$ pairs/cm
(`ion_track_generation_2d`, `etna/core/single_particle.py`; `track_pairs`,
`etna/core/track/phd.py`); `ion_let_profile` gives $L(z)$ along a Bragg path from the same
NIST tables. C-12 in the MC sampler uses $z_{\rm eff}^2$ times PSTAR at the same velocity (15.1).

### 8.5 Track structure (radial distribution of the pairs)

The 2D single-particle tables (10.3) use a vertical sheet with a Gaussian lateral profile of
$\sigma = 1$ µm. The track engine (`etna/core/track/radial.py`) uses the cumulative fraction
$F(r)$ of pairs within radius $r$ and lumps the pairs of each annulus into its finite volume
exactly, $F(r_b)-F(r_a)$, so the total is conserved on any mesh:

- `gaussian(w)`: $F = 1-e^{-r^2/w^2}$ (the generic commercial-TCAD heavy-ion form; Sentaurus'
  exact default is not verified);
- `core_penumbra(r_min, r_max)`: constant density inside $r_{\min}$, $\propto1/r^2$ to $r_{\max}$
  (the Local Effect Model form, Korcyl arXiv:1410.5250 eq. 1.20 [KOR14]);
- `whk(beta)`: the Waligórski-Hamm-Katz radial dose with water constants from Korcyl eqs.
  2.3-2.12, transferred to SiC by mass-thickness scaling ($\times\rho_{\rm w}/\rho_{\rm SiC}$),
  **our approximation**, normalised to the stopping-power LET;
- `thermalized(model, s)`: the model convolved with a planar Gaussian of width $s$. The
  thermalisation length of 4H-SiC is **not known**, so every result that depends on the core
  density is conditional on $s$.

**Status:** NIST data and their Bragg composition: **literature data + assumption** (additivity);
the generation formulas are exact bookkeeping (**verified** by the worked example); the radial
models are **assumptions** with literature shapes.

**Worked example** (`worked_examples.py generation`):

| quantity | code | by hand |
| --- | --- | --- |
| $G$ per 1 Gy/s in SiC | $2.5689\times10^{15}$ cm$^{-3}$ s$^{-1}$ | $3.21\times10^4/(7.8\times1.602\times10^{-12})$ = $2.5689\times10^{15}$ |
| $L_{\rm SiC}$, proton 5.6 / 62 / 150 MeV | 18.2136 / 2.7761 / 1.4511 keV/µm | $(0.7004S_{\rm Si}+0.2996S_{\rm C})\rho$: identical |
| pairs per µm, same | 2335.1 / 355.9 / 186.0 | $L/E_{\rm pair}$ |
| $G$, 62 MeV, $10^8$ p/cm$^2$/s (fluence mode) | $3.5591\times10^{14}$ cm$^{-3}$ s$^{-1}$ | $10^8\times277.61\times10^6/7.8$ |
| 5.486 MeV alpha: CSDA range, Bohr $\sigma_R$ | 18.429 µm, 0.165 µm | |
| alpha peak / entrance, peak depth | 2.488 at 16.94 µm | |
| alpha pairs | 703 333 | $5.486\times10^6/7.8$ |

## 9. Signal formation: conduction, displacement and induced current

### 9.1 Terminal current

The current at a contact from the continuity equations is the **conduction** current
(`extract_contact_current`, `etna/core/drift_diffusion.py`). In a transient the terminal current
also contains the **displacement** current, the time derivative of the contact charge of the
Poisson equation, $Q_c = \oint\varepsilon\mathbf E\cdot d\mathbf S$ (`extract_contact_charge`).
For an implicit Euler step of length $\Delta t_k$ (`displacement_current`)

$$
I_{\rm disp,k} = \frac{Q_{c,k}-Q_{c,k-1}}{\Delta t_k},\qquad I_{\rm tot} = I_{\rm cond} + I_{\rm disp}.
$$

$\sum_kI_{\rm disp,k}\Delta t_k = Q_{c,\rm end}-Q_{c,\rm start}$ exactly, so conduction and total
current carry the same charge once the device returns to its initial state: the collected charge
and CCE are right with either; the pulse **shape** needs the total. In 1D the totals at the two
contacts sum to zero at every step (Kirchhoff); the conduction parts alone do not.

### 9.2 Shockley-Ramo in 1D

With contacts at both ends of a 1D device of length $L_d$ at fixed biases, the weighting field
of either contact is $1/L_d$ everywhere, so the induced current is the device-averaged
conduction current (`ramo_current_1d`):

$$
I_{\rm ind} = \frac{1}{L_d}\int_0^{L_d}(J_n+J_p)\,dx ,
$$

an independent route to $I_{\rm cond}+I_{\rm disp}$ (it uses all edge currents, not the contact
equations). Both are checked against each other and against the parallel-plate result
$I = qNv/d$ (`tests/test_displacement_current.py`): **verified**. 2D weighting fields are not
computed; the 2D pulses use the terminal currents. The Allpix Squared export writes the
planar-pad weighting potential $\phi_w = z/t$ (`export_allpix_fields`, `etna/core/mc_exchange.py`:
exact for a pad much wider than thick; the files are not checked against an Allpix run).

### 9.3 Time integration of the devsim transients

`simulate_single_particle` (`etna/core/single_particle.py`) injects the pairs as a generation
pulse $G/\Delta t_{\rm inj}$ during one BDF1 step ($\Delta t_{\rm inj} = 1$ ps), then steps
$\Delta t = \mathrm{clip}(0.1\,t,\ 10^{-13},\ 10^{-9})$ s with BDF1 until the charge added per step
falls below $10^{-7}$ of the total for three steps. The collected charge is the backward-Euler
sum, exact for the time-discrete system (`collected_charge`):

$$
Q_{\rm col} = \sum_k\big(|I_k|-|I_{\rm dark}|\big)\,\Delta t_k .
$$

The FLASH pulse transient (`TransientSolver`, `etna/core/transient.py`) is BDF1 on a trapezoidal
pulse envelope with phase-based steps (rise/fall $t_{\rm rise}/10$, plateau $\Delta t_{\max}$).
**Deviation:** devsim's charge-based step rejection is disabled there (`charge_error=1e10`), so
neither devsim transient has a local error control; their accuracy rests on the step rules. The
track engine has one (TR-BDF2, 16.3).

## 10. Charge collection efficiency

### 10.1 Definitions

ETNA uses several CCEs; they are not interchangeable.

| name | definition | code | used by |
| --- | --- | --- | --- |
| steady-state DD CCE | $\mathrm{CCE} = \lvert I_{\rm col}\rvert/I_{\rm gen}$, $I_{\rm gen} = q\int G\,dx$ (trapezoid = box rule in 1D), $I_{\rm col}$ the steady contact current in the balanced node-set form | `compute_cce_from_dd`, `etna/core/charge_collection.py`; engine `etna/core/batched/sweeps.py` | `run_cce`, CCE vs fluence, temperature |
| 2D steady CCE | the same with $I_{\rm gen}$ the box-rule integral over **non-contact** nodes (`injected_generation`) | `compute_cce_2d`, `etna/core/charge_collection_2d.py` | lateral scans, microdosimetric sweep |
| single-particle CCE | $Q_{\rm col}/Q_{\rm gen}$ of a transient (9.3) | `simulate_single_particle` | CCE(LET) tables |
| track CCE | $1-(\text{recombined}+\text{wrong-contact minority}+\text{stored})/N_{\rm gen}$ from the pair ledger | `etna/core/track/solver.py` `TrackResult.cce` | track engine |
| Hecht | closed form, 6.4 | `hecht_cce` | analytic reference |

The steady-state CCE is the ratio of two currents under a constant, low-injection generation
(peak $10^{18}$ cm$^{-3}$ s$^{-1}$, so the excess density in the depleted epi is ~$10^8$
cm$^{-3}$, far below $N_D$). It measures the loss by recombination and by diffusion to the wrong
contact; it is **not** a single-particle CCE: SRH in the depletion region needs both carriers at
once, which a steady background supplies and a single pulse does not (USAGE §19). The direct
majority-carrier current at an n$^+$ contact is a small difference of large SG terms dominated by
round-off, so the steady current is read in a balanced form (`extract_steady_state_current`: the
volume generation minus the fluxes leaving through the other contacts); the transients integrate
the anode current when there is an n$^+$ layer.

**Deviation, fixed and documented.** Before 2026-10-18 the 2D $I_{\rm gen}$ counted generation on
contact nodes that the solver never injects, a mesh-dependent CCE deficit (1.26e-3 on the legacy
planar diode at the cathode).

### 10.2 What limits the CCE in the default device

Generation in the epi only; holes generated in the undepleted epi diffuse; the n$^-$/n$^+$
high-low junction reflects holes; the p$^+$ is a dead layer only if generation is put there.
Midgap SRH with $\tau_n = 1$ ns, $\tau_p = 600$ ns. With these inputs the CCE of the default
device is $>0.999$ above 10 V: a **prediction**, not checked against a measured CCE of this
device. The microdosimeter's measured alpha CCE-vs-bias shape is **not** reproduced with the
library lifetimes (USAGE §19.3).

### 10.3 Single-particle CCE(LET) tables

`build_cce_let_table` (`etna/core/single_particle.py`) runs the 2D transient of 9.3 for a vertical
sheet track through the epi (Gaussian lateral $\sigma = 1$ µm, centred on the symmetry plane of the
half device) over a log grid of LET and stores $\mathrm{CCE}(L)$, $Q_{\rm gen}$, $Q_{\rm col}$
and the collection time; `method="linear"` synthesises low LETs from the highest one once a
simulated anchor agrees to $2.5\times10^{-7}$. The shipped tables `data/cce_let_table_100um.json`
and `_300um.json` (10 µm epi, 50 V, `c_axis`) hold 20 LETs from 0.5 to 500 keV/µm. The
Monte Carlo coupling applies them event by event (`process_mc_ensemble`,
`etna/core/mc_response.py` `pulse_heights`: $Q = \mathrm{CCE}(L_k)\,\epsilon_k/E_{\rm pair}$, LET
per event from its chord or the thickness).

**Limit.** The 2D sheet track at these densities stays linear: the table is CCE at low injection
and cannot show plasma recombination; the pulse-height defect needs the track engine (10.4).

**Worked example** (`worked_examples.py cce`): `run_cce(DeviceConfig())` gives 0.999433 /
0.999706 / 0.999799 / 0.999847 at $-10$ / $-20$ / $-30$ / $-40$ V (batched engine). The 100 µm
table spans CCE 0.999966-0.999974 over 0.5-500 keV/µm, and $Q_{\rm gen}/L = 0.10385$ fC per
keV/µm at every LET (constant to $7\times10^{-17}$, i.e. linear generation): by hand, half of the
track (the mirror plane cuts the Gaussian) over 10 µm is
$\tfrac12\times10\,\text{µm}\times1\,\text{keV/µm}/7.8\,\text{eV}\times q = 0.1027$ fC; the
1.1 % excess is the box-rule volume of the boundary nodes at the junction and the epi end
(the stored $Q$ is per unit length along the sheet, labelled fC).

### 10.4 Single-ion tracks, plasma recombination and the pulse-height defect (opt-in)

The track engine (`etna/core/track/`, facades in `etna/api/track.py`) solves one ion at normal
incidence in $(r,z)$ with the radial models of 8.5, SRH + Auger (constant or injection-dependent),
optional carrier-carrier scattering (Conwell-Weisskopf with the Lades 4H-SiC fit, its derivative
left out of the Jacobian), instantaneous or Gaussian-pulsed deposition, and the dark state of the
batched 1D engine. Its pair ledger (generated = recombined + stored + collected) closes to the
Newton tolerance; the pulse-height defect it reports is the **recombination part only** (no
window, dead-layer-by-geometry or nuclear-stopping defect).

A tilted track is composed (`etna/core/track/tilted.py`) as

$$
\mathrm{CCE}(V) = \mathrm{CCE}_{\rm column}(V) - \sum_kw_k\,f_{\rm ex}\big(s_k,E_\perp(z_k;V)\big),\qquad f_{\rm ex} = f(L(s),E_\perp)-f(L\to0,E_\perp),
$$

with $\mathrm{CCE}_{\rm column}$ the linear response of the tilted generation on the 1D column,
and $f$ the recombined fraction of a long, locally uniform track in the transverse field
$E_\perp = |E|\sin\theta$ (`cross_section.slice_loss`, 2D Cartesian half plane, absorbing far
edges). Stated approximations: $E\cos\theta$ dropped (overestimates the loss), box half-width $R$
matters at low field, undepleted segments use the $E\to0$ depleted excess
([USAGE §19.2](USAGE.md#192-tilted-tracks-transverse-slices-and-the-25d-composition)).

**Status:** numerics **verified** (ledger, energy theorem, low-injection limit equals the 1D
engine, `tests/test_track_engine.py`, `tests/test_track_tilted.py`); PHD and plasma loss
**predicted**, conditional on the core radius; the tilted fit of the Fig. 4b data is a
**calibration** of the lifetimes (USAGE §19.3), not a reproduction.

**Worked example** (`worked_examples.py tracks`): `gaussian(50 nm)` puts $F$ = 4.0e-4 / 0.039 /
0.982 / 1.000 of the pairs within 1 nm / 10 nm / 100 nm / 1 µm; one track of LET 1000 keV/µm at
50 V with that core gives CCE 0.999872, a recombined fraction $1.28\times10^{-4}$ and a ledger
residual $1.1\times10^{-10}$ (`MESH_FAST`).

## 11. FLASH and high injection

`cce_vs_dose_rate` (`etna/core/flash_recombination.py`; engine `etna/core/batched/dark.py`
`generation_cce`) solves the steady state with midgap SRH + Auger (6.2) under the flat proton
generation of 8.4 for each dose rate, and $\mathrm{CCE} = |J|/(q\int G\,dx)$.
`solve_with_continuation` ramps the generation for convergence on the devsim path.

**What it can and cannot say.** At 20-230 Gy/s ($G\approx5\times10^{16}$-$6\times10^{17}$
cm$^{-3}$ s$^{-1}$) the excess densities are far too low for high-injection Auger to matter in the
epi; the only Auger effect is a **dose-rate-independent** low-injection offset from
minority-hole recombination in the $10^{19}$ n$^+$, where $1/(C_nN_D^2)\approx20$ ns $<\tau_p$.
The genuine high-injection physics (e-h plasma screening of the field, ambipolar transport,
conductivity modulation, track funnelling) is **not implemented** in this path; radiative
recombination is not wired. **Status:** a **numerical sensitivity bound**, not a FLASH
prediction (README); the track engine (10.4) is the place where single-ion plasma effects are
modelled.

**Worked example** (`worked_examples.py flash`): CCE(20 Gy/s) = 0.999810726, CCE(230 Gy/s) =
0.999810740 (a change of $1.4\times10^{-8}$), and the SRH-only reference 0.999811068: the Auger
offset is $3.4\times10^{-7}$, independent of dose rate.

## 12. Radiation damage

The damage chain is: particle fluence $\Phi$ at energy $E$ $\to$ equivalent fluence
$\Phi_{\rm eq} = \kappa(E)\,\Phi$ (12.1-12.2) $\to$ defect densities $N_k = \eta_k\Phi_{\rm eq}$
(12.3) $\to$ either a scalar donor removal (legacy, 12.4) or charged deep traps in Poisson with
steady or kinetic occupancy plus a calibrated direct removal (default, 12.5) $\to$ trapping
lifetimes (12.6) $\to$ C-V, CCE and leakage at fluence. Entry point:
`run_radiation_damage` (`etna/api/simulation.py`), model resolution `damage_model_of` /
`resolve_damage_model` (`etna/core/charge_collection.py`).

### 12.1 NIEL and the 1 MeV neutron equivalent

The displacement damage of a particle is assumed proportional to its non-ionising energy loss
(NIEL, MeV cm$^2$/g), and fluences of different particles and energies are compared through the
hardness factor

$$
\kappa(E) = \frac{\mathrm{NIEL}_x(E)}{\mathrm{NIEL}_{\rm ref}},\qquad \Phi_{\rm eq} = \kappa(E)\,\Phi .
$$

Data (`etna/core/niel_tables.py`, provenance in `data/niel/README.md`):

- proton NIEL of Si and of SiC (compound, Si:C = 1:1, $E_d$ = 21/21 eV) from the SR-NIEL web
  calculator, version SR-NIEL-7 ver. 11.1, Hadron + Coulomb, no hadronic scaling, Robinson
  partition, fetched 2026-10-20 [SRNIEL]; 1-1000 MeV used, log-log interpolation, no
  extrapolation (`get_hardness_factor` raises outside);
- the ASTM E722 normalisation of silicon, $D_n(1\ \text{MeV}) = 95$ MeV mb [FRE02]:

$$
\mathrm{NIEL}_{n,\rm Si}(1\ \text{MeV}) = \frac{95\times10^{-27}\,\text{MeV cm}^2\times N_A}{28.0855\ \text{g}} = 2.0370\times10^{-3}\ \text{MeV cm}^2/\text{g};
$$

- neutron NIEL of Si and 4H-SiC computed with NJOY2016 HEATR on ENDF/B-VIII.0 (Lindhard /
  Robinson partition, $E_d$ 21 eV) by `scripts/build_neutron_niel.py`, stored in
  `data/niel/neutron_*.csv(.gz)` and read by `etna/core/neutron_niel.py`:
  $\mathrm{NIEL}_{n,\rm SiC}(1\ \text{MeV}) = 2.29889\times10^{-3}$, Watt-fission-spectrum
  averages $2.18058\times10^{-3}$ (Si) and $2.33030\times10^{-3}$ (SiC), hence

$$
r_n = \frac{\langle\mathrm{NIEL}_{\rm SiC}\rangle_{\rm fission}}{\langle\mathrm{NIEL}_{\rm Si}\rangle_{\rm fission}} = 1.0687 .
$$

  Checks recorded with the data: Si pointwise at 1 MeV = 0.998 of the ASTM value; the Si fission
  average within 0.7 % of the Griffin damage function; an independent elastic integration
  reproduces HEATR to $10^{-4}$; ENDF/B-VII.1 changes $r_n$ by $-0.16$ %; reactor-like spectra
  give $r_n$ 1.02-1.20. No published SiC 1 MeV value was found to compare with.

### 12.2 The three NIEL conventions

The introduction rates of the trap table are quoted per 1 MeV neutron equivalent *as labelled
by reactor-neutron data* (Burin et al. [BUR24t] fit TRIGA data labelled in neq [GSP23]), and the
same group converts protons with **silicon** NIEL (Radmanovac et al. [RAD25]: $k\approx0.9$ at
252.7 MeV; SR-NIEL gives 0.855). ETNA offers three conventions
(`hardness_table(particle, convention)`, `etna/core/radiation_damage.py`):

$$
\kappa_{\text{si-neq}}(E) = \frac{\mathrm{NIEL}_{x,\rm Si}(E)}{\mathrm{NIEL}_{n,\rm Si}(1\,\text{MeV})},\qquad
\kappa_{\text{sic-damage}}(E) = \frac{\mathrm{NIEL}_{x,\rm SiC}(E)}{r_n\,\mathrm{NIEL}_{n,\rm Si}(1\,\text{MeV})},\qquad
\kappa_{\text{sic-1mev}}(E) = \frac{\mathrm{NIEL}_{x,\rm SiC}(E)}{\mathrm{NIEL}_{n,\rm SiC}(1\,\text{MeV})}.
$$

**Derivation of the default ("sic-damage", since 2026-10-21).** Suppose the introduction rates
$\eta$ were measured on SiC with a reactor (fission-like) neutron spectrum whose fluence was
labelled in Si-neq: $\Phi_{\rm eq}^{\rm label} = \Phi_n\langle\mathrm{NIEL}_{\rm Si}\rangle/\mathrm{NIEL}_{n,\rm Si}(1)$.
The SiC damage those neutrons did is $\propto\Phi_n\langle\mathrm{NIEL}_{\rm SiC}\rangle = r_n\,\mathrm{NIEL}_{n,\rm Si}(1)\,\Phi_{\rm eq}^{\rm label}$.
A proton fluence $\Phi_p$ does SiC damage $\propto\Phi_p\,\mathrm{NIEL}_{p,\rm SiC}(E)$. Equating the
damage gives the Si-labelled equivalent fluence that produces the *same SiC damage*,

$$
\Phi_{\rm eq} = \frac{\mathrm{NIEL}_{p,\rm SiC}(E)}{r_n\,\mathrm{NIEL}_{n,\rm Si}(1)}\,\Phi_p = \kappa_{\text{si-neq}}(E)\,\frac{r_p(E)}{r_n}\,\Phi_p,\qquad r_p = \frac{\mathrm{NIEL}_{p,\rm SiC}}{\mathrm{NIEL}_{p,\rm Si}} ,
$$

so the per-neq rates carry over unchanged (a reactor-neutron fluence has the same $\kappa$ in
"si-neq" and "sic-damage") and the SiC damage per proton drops by $r_p/r_n$ ($\times0.80$ at
62 MeV). "si-neq" assumes $r_p = r_n$ and overstates the SiC damage per proton by $r_n/r_p$.
"sic-1mev" is SiC-native, for rates quoted per SiC 1 MeV neq (none of the shipped ones are).
Proton-calibrated numbers (12.5) are per-proton facts and are re-fitted per convention.

**Status:** NIEL tables **literature data**; the SiC neutron normalisation is **computed, not
checked against a published SiC value**, and the issue stays open (#8); the convention is a
maintainer decision. Known systematic: measured Si proton hardness factors are 35-50 % below the
NIEL prediction at 10-25 MeV (arXiv:2408.12410, cited in `radiation_damage.py`), i.e. at and
below the 5.6 MeV API default.

**Worked example** (`worked_examples.py niel`): at 62 MeV, $\mathrm{NIEL}_{p,\rm Si} = 3.95580\times10^{-3}$
 and $\mathrm{NIEL}_{p,\rm SiC} = 3.37920\times10^{-3}$ MeV cm$^2$/g give
$\kappa$ = 1.9420 (si-neq), 1.5523 (sic-damage) and 1.4699 (sic-1mev), equal to the hand ratios
above; $10^{12}$ p/cm$^2$ is $1.942$, $1.552$ and $1.470\times10^{12}$ neq/cm$^2$.

### 12.3 Defects and their introduction

$N_k = \eta_k\Phi_{\rm eq}$ (`defect_concentration`). The trap table `RadiationDamageParams`
(Burin et al. [BUR24t] Table I; the rates cited there to Hazdra et al., not opened):

| level | type | $E_c-E_t$ (eV) | $\eta$ (cm$^{-1}$) | $\sigma_n$ (cm$^2$) | $\sigma_p$ (cm$^2$) |
| --- | --- | --- | --- | --- | --- |
| Z1/2 | acceptor | 0.67 | 5.0 | $2\times10^{-14}$ | $3.5\times10^{-14}$ |
| EH4 | acceptor | 1.03 | 2.4 | $5\times10^{-13}$ | $5\times10^{-14}$ |
| EH6/7 | donor | 1.60 | 1.6 | $1.4\times10^{-14}$ (see note) | $3.8\times10^{-14}$ |

**Deviation (EH6/7 $\sigma_n$).** Burin's Table I prints $9\times10^{-12}$; ETNA keeps the audit
value $1.4\times10^{-14}$ in this table and uses the emission-consistent $10^{-13}$ in the
kinetic presets. Neither has a primary *capture* measurement: with $E = 1.60$ eV, $9\times10^{-12}$
emits 40-160$\times$ faster and $1.4\times10^{-14}$ 4-15$\times$ slower than the measured EH6/7 DLTS
rates (Hazdra et al. 2016, Danno & Kimoto 2006, Alfieri & Kimoto 2013, as opened in
`lab/radiation_traps/notes/MECHANISMS.md`).

Annealing (opt-in, `annealing_fraction`): first-order Arrhenius recovery
$f = 1-\exp(-\nu_0te^{-E_a/k_BT})$ with $\nu_0 = 10^{13}$ Hz and $E_a$ = 4.5 (Z1/2), 1.8 (EH4),
3.2 eV (EH6/7), **estimated** from thermal-stability onsets, not measured Arrhenius analyses.

### 12.4 Legacy scalar carrier removal (`radiation_traps=False` / `"legacy"`)

$$
N_D(x;\Phi) = \max\big(N_D(x)-\eta_{\rm rem}\Phi_{\rm eq},\ 0\big),\qquad \eta_{\rm rem} = 5.0\ \text{cm}^{-1},
$$

on the epi nodes (`apply_carrier_removal`). Full compensation is reached at
$\Phi_{\rm comp} = \min_xN_D/(\eta_{\rm rem}\kappa)$; the doping-based model is ill-posed beyond it,
so those points are returned as NaN (`full_compensation_reason`). Measured removal at 252.7 MeV
is 5-7.5 per neq [RAD25].

### 12.5 Charged deep traps in Poisson (default: calibrated kinetic preset)

`etna/core/radiation_traps.py` replaces the removal scalar by the levels themselves as charged
species. On the epi nodes $N_k = \eta_k\Phi_{\rm eq}$ and the trap charge enters Poisson (3.1):

$$
\rho_{\rm trap}/q = \sum_kN_k\,(d_k - f_k),\qquad d_k = 1\ \text{(donor)},\ 0\ \text{(acceptor)} ,
$$

(an acceptor carries $-N_kf_k$, a donor $+N_k(1-f_k)$; `trap_charge`). The donors $N_D(x)$ are not
reduced by the traps themselves: removal *emerges* from the electrons the acceptors bind.

**Rates.** $c_{n,p} = \sigma_{n,p}v_{\rm th}$, $v_{\rm th} = \sqrt{3k_BT/m^*}$ ($m^*$ = 0.77, 1.0
$m_0$), $e_n = c_nn_1$, $e_p = c_pp_1$,

$$
n_1 = N_C(T)\,e^{-(E_c-E_t)/k_BT},\qquad p_1 = n_i^2/n_1 ,
$$

with the unanchored $N_C$ and the solver's anchored $n_i$ (so that at equilibrium $f$ is exactly the
Fermi-Dirac occupancy $n/(n+n_1)$; this makes $e_p$ about 0.42 of an $N_V$-based value). Degeneracy
factor 1 (**assumption**: Table I gives none).

**Steady occupancy** (`radiation_traps=True`, the Burin table):

$$
f = \frac{c_nn + e_p}{c_n(n+n_1) + c_p(p+p_1)} ,
$$

and for a two-level negative-U trap (Z1, Z2 of the kinetic presets) the zero-net-flux chain
$0\leftrightarrow1\leftrightarrow2$ electrons (`occupancy`). The dark state at the applied bias is
solved self-consistently with Poisson (Newton on $\psi$ with the trap charge and its derivative,
step clamp 2 $V_t$); for a CCE point the dark trap charge is then **frozen** in the illuminated
solve (a real alpha pulse does not repopulate traps with ms to $10^{14}$ s time constants).
Known consequence: in depletion the steady state empties both acceptors and leaves EH6/7 ~90 %
ionized, so the depleted $N_{\rm eff}$ *grows* with fluence, the wrong sign against measured C-V
removal; the model never type-inverts.

**Kinetic occupancy** (`etna/core/trap_kinetics.py`): the occupancy integrates the history,

$$
\frac{df}{dt} = c_nn(1-f) - e_nf - c_ppf + e_p(1-f) = a - \lambda f,\quad a = c_nn+e_p,\ \lambda = c_n(n+n_1)+c_p(p+p_1),
$$

with backward Euler folded into the Poisson Newton, $f_{\rm new} = (f_{\rm old}+\Delta t\,a)/(1+\Delta t\,\lambda)$
($\Delta t\to\infty$ is the steady kernel, $\Delta t = 0$ the frozen charge), quasi-static carriers
(dark state with flat quasi-Fermi levels, optional photo-carriers), and adaptive $\Delta t$ by step
doubling on $\max|f(\Delta t)-f(\Delta t/2,\Delta t/2)|$. C-V observables freeze the trap charge at
the AC frequency; $C(f)$ of a semi-insulating epi uses the dielectric-relaxation slab
$Y = 1/\sum_eh_e/(\sigma_e+i\omega\varepsilon)$, $\sigma = q(\mu_nn+\mu_pp)$. Not modelled: the trap
response at the AC frequency, transport-limited quasi-Fermi splitting in a semi-insulating epi,
Poole-Frenkel enhancement, temperature changes within a run.

**The default preset** `radmanovac_calibrated` (`etna/core/trap_presets.py`, selected by
`radiation_traps="auto"` inside the batched engine's scope; outside it "auto" falls back to the
legacy removal with a logged note):

- `literature_kinetic` (**predicted**, nothing fitted): Z1/Z2 as two-level negative-U acceptors
  with Capan et al. 2018 Table I energies and Arrhenius emission (Z1: 0.48 / 0.59 eV,
  Z2: 0.41 / 0.67 eV, $[{\rm Z2}]:[{\rm Z1}] = 4.4$) [CAP18], 2.5 carbon vacancies per neq
  (**assumption**: the DLTS concentration counts two emitted electrons per vacancy), EH4 from
  Table I, EH6/7 with $\sigma_n = 10^{-13}$; kinetic history: 0 V equilibrium start, 1 V/s ramp,
  60 s hold at bias for CCE; occupancy-weighted lifetimes (12.6). The emission cross-section is
  converted from Capan's Arrhenius prefactor so that the 300 K rates are Capan's
  (`arrhenius_sigma`: $e_n = AT^2e^{-E_a/k_BT} = \sigma v_{\rm th}N_Ce^{-E_a/k_BT}$);
- plus a **direct donor removal** $N_D\to\max(N_D-d\,\Phi_{\rm eq},0)$, **calibrated** to the three
  C-V slopes of Radmanovac et al. Table 2 (252.7 MeV protons, 1 MHz C-V; 5.60 ± 0.81, 4.60 ± 0.33,
  4.48 ± 0.26 per p/cm$^2$ on onsemi and CNM devices, *not* this detector) [RAD25]. Because the
  slopes are per proton, $g = \kappa_c(252.7)\,[G_{\rm lit}(\kappa_c)+d_c]$ is re-fitted in each
  convention $c$: $d$ = 2.87 ± 0.22 (si-neq), **4.69 ± 0.30 (sic-damage)**, 5.10 ± 0.32 (sic-1mev)
  per neq ($\chi^2$ about 1 for 2 dof; the ± are data spreads over depletion widths). A slow deep
  acceptor near midgap would give the same room-temperature C-V; `deep_acceptor=True` carries the
  same rate that way. **The removal mechanism is not identified.**

In the neutral epi every acceptor is filled, so the removal rate is
$\eta_{\rm neutral} = 2\times2.5\ (\text{Z1+Z2, two electrons each}) + 2.4\ (\text{EH4}) + d$ per neq
(`neutral_removal_rate`), and the compensation fluence is
$\Phi_{\rm comp} = \min_xN_D/(\eta_{\rm neutral}\,\kappa)$; beyond it the semi-insulating epi is
still solved.

**Status:** levels and cross-sections **literature values** (with the stated EH6/7 exception);
removal **calibrated** on other devices at one energy and one temperature; CCE and C-V vs fluence
of this device **predicted** ([USAGE §6](USAGE.md#6-the-run_-facades)). `trap_calibration.py`
fits multiplicative scales on the unmeasured introduction rates to measured C-V curves of *this*
device when they exist (Gauss-Newton covariance in log space, identifiability flags), with the
geometry and pristine $N_D$ fixed from measurement.

### 12.6 Trapping lifetimes

Legacy (`compute_K_tau`, `degraded_lifetime`):

$$
\frac{1}{\tau(\Phi)} = \frac{1}{\tau_0} + K_\tau\Phi_{\rm eq},\qquad K_{\tau,n} = v_{\rm th,n}\sum_k\eta_k\sigma_{n,k},\quad K_{\tau,p} = v_{\rm th,p}\sum_k\eta_k\sigma_{p,k},
$$

i.e. every level counted as fully available for capture (an optional "logarithmic" form
$\tau = \tau_0/(1+K_\tau\tau_0\Phi)^{0.8}$ is empirical). With the trap presets
(`lifetimes="occupancy"`, `occupancy_lifetimes`) only the empty share captures electrons and only
the filled share captures holes:
$1/\tau_n = 1/\tau_{0,n} + \langle\sum_kN_k\sigma_{n,k}v_{\rm th,n}(1-f_k)\rangle_{\rm epi}$ (and the
mirror for holes, with the multi-level availabilities of `capture_availability`), averaged over the
epi nodes. The lifetimes stay region-wide scalars (midgap SRH, 6.1).

**Worked example** (`worked_examples.py damage`): $v_{\rm th,n} = 1.3310\times10^7$ cm/s and
$K_{\tau,n} = (5.0\cdot2\times10^{-14} + 1.6\cdot1.4\times10^{-14} + 2.4\cdot5\times10^{-13})\,v_{\rm th} = 1.7601\times10^{-5}$
 cm$^2$/s (code = hand), $K_{\tau,p} = 4.1555\times10^{-6}$ cm$^2$/s. Emission
at 300 K, $e_n = \sigma_nv_{\rm th}N_Ce^{-E/k_BT}$: Z1/2 25.1 /s (1/e = 40 ms), EH4
$5.6\times10^{-4}$ /s (30 min), EH6/7 $4.2\times10^{-15}$ /s ($2.4\times10^{14}$ s). Compensation of the
default stack by 62 MeV protons: $1.802\times10^{12}$ p/cm$^2$ with the default model
($3.3825\times10^{13}/(12.09\times1.5523)$, $12.09 = 5+2.4+4.69$), $3.484\times10^{12}$ with the
legacy removal in si-neq ($3.3825\times10^{13}/(5\times1.942)$).

## 13. Temperature and self-heating

### 13.1 Temperature as a parameter

Every temperature dependence in one place (the code sites are in chapters 2, 5, 6 and 12):
$E_g(T)$ Varshni; anchored $n_i(T)$; $N_C,N_V\propto T^{3/2}$; $V_t = k_BT/q$;
$\tau\propto T^{1.72}$; the mobility laws of each set; $N_t\propto n_i(T)$ in the dark current;
the legacy Al fraction $\times(1+0.5(T-300)/300)$; the Boltzmann factors of the opt-in ionization
and of the trap emission. Not temperature-dependent: $\varepsilon_r$, the ionization $\Delta E_0$,
$\alpha$, $g$, the Auger coefficients, $\kappa_{\rm th}$, $c_p$. `DeviceConfig` accepts 100-600 K;
the C-V engine is checked at 250-400 K and the dark-current engine at 150-450 K
(`etna/core/temperature_sweep.py`, `run_temperature_sweep`). **Status:** the shapes are
**predicted**, assumption-driven ("bias and temperature shapes are model assumptions", README).

### 13.2 Lattice heating of a single ion track (opt-in)

`etna/core/track/thermal.py` (`LatticeHeat`), switched on with `thermal=` in the track engine:

$$
C_v\frac{\partial T}{\partial t} = \nabla\cdot(\kappa_{\rm th}\nabla T) + H,\qquad
H = \mathbf J\cdot\mathbf E + E_g\,(U-U_{\rm dark}) + \text{prompt heat},
$$

with $\mathbf J\cdot\mathbf E$ the SG edge power $I_e(\psi_k-\psi_l)$ split half to each node (the
discrete form whose integral closes the electrical energy theorem), $E_g$ per excess recombination,
and the prompt heat $(E_{\rm pair}-E_g)$ per deposited pair added as an instantaneous temperature
rise. Constants held at their 300 K values: $\kappa_r = 4.71$ W/(cm K) (basal plane, radial),
$\kappa_z = 3.24$ W/(cm K) (along c), semi-insulating 4H-SiC, TDTR, Qian, Jiang & Yang
[QIA17] Table 1; $c_p = 0.69$ J/(g K) (Ioffe, the 6H value used for all polytypes; opened for
this value [IOFFE]), $\rho = 3.21$ g/cm$^3$, so $C_v = 2.2149$ J/(cm$^3$ K). Boundaries adiabatic
(default) or the cathode plane held at $T_0$. The heat equation uses the electrical step's TR-BDF2
stages and weights, so the energy ledger closes exactly. With `couple=True` mobility and SRH
lifetimes follow the node temperature after each step (first-order splitting).

**Deviation from a full electrothermal model:** transport stays isothermal SG (no thermodiffusion,
no Seebeck current), $\kappa_{\rm th}(T)$ and $c_p(T)$ are constant, and no other solver has a heat
equation (the dosimeter's steady self-heating is not modelled; at 18 pA and 60 V it would be
~1 nW).

**Status:** numerics **verified** (analytic line source within 3 %, ledger $\le10^{-7}$,
`tests/test_track_engine.py`); constants **literature values**; results **predicted** (USAGE
§19: 1.8 K for a 50 nm core, 118 K for picoseconds and $\Delta$CCE $-2.5\times10^{-4}$ for a 10 nm
core at LET 1000 keV/µm).

**Worked example** (`worked_examples.py device`, section 7): the analytic line source
$\Delta T(0,t) = E_L/\big(\pi C_v(w^2+4at)\big)$, $a = \kappa/C_v$, for the whole deposit of
LET 1000 keV/µm with $w$ = 10 nm, gives 24.22 / 2.68 / 0.27 K at 1 / 10 / 100 ps in the code and
by hand (the engine's prompt heat is only the $(E_{\rm pair}-E_g)/E_{\rm pair}$ share of it).

## 14. Microdosimetry

### 14.1 Lineal energy and its distributions

Following the ICRU Report 36 conventions as restated in `etna/core/microdosimetry.py` (the report
itself was not opened by this project), the lineal energy of an event depositing $\epsilon$ in the
sensitive volume (SV) is $y = \epsilon/\bar l$ (keV/µm). On log-spaced bins (`make_y_bins`,
default $10^{-2}$-$10^4$ keV/µm, 50 per decade, centres the geometric means)
`lineal_energy_spectrum` computes

$$
f(y_j) = \frac{N_j}{N\,\Delta y_j},\qquad y_F = \sum_jy_jf(y_j)\Delta y_j,\qquad d(y_j) = \frac{y_jf(y_j)}{y_F},\qquad y_D = \sum_jy_jd(y_j)\Delta y_j ,
$$

and plots $y\,f(y)$ and $y\,d(y)$ on a log axis. **Deviation:** these $y_F$, $y_D$ are binned
estimates (bin centres), and $N$ counts every event, including those outside the bin range, so
$\int f\,dy<1$ when events fall outside (a warning is logged above 1 %). The event-level
estimators of `etna/core/mc_response.py` `lineal_moments` are exact and come with standard
errors:

$$
y_F = \frac{\langle\epsilon\rangle}{\bar l},\qquad y_D = \frac{\langle\epsilon^2\rangle}{\langle\epsilon\rangle\,\bar l},\qquad
\mathrm{var}(y_D) \approx \frac1N\left[\frac{v_2}{m_1^2}+\frac{m_2^2}{m_1^4}v_1-\frac{2m_2}{m_1^3}c_{12}\right]
$$

(delta method, $m_{1,2}$ the sample means of $y$, $y^2$, $v$ their variances, $c_{12}$ their
covariance); `adaptive_lineal_energy` samples until both relative errors are below a target.

### 14.2 Mean chord length

`mean_chord_length`: Cauchy's theorem for a convex body under isotropic (µ-random) irradiation,
$\bar l = 4V/S$, for a rectangular SV $w\times t\times d$:
$\bar l = 4wtd/[2(wt+wd+td)]$; without width and depth the slab limit $\bar l = 2t$.
`run_microdosimetry(chord_model="cauchy" | "slab")` chooses (depth defaults to the width).
For a parallel beam the event chord is $t/\cos\theta$ (`mc_response`, `irradiation="normal"` or an
angle), and the MC sampler returns each event's own chord. Using the isotropic $\bar l$ for a
normally incident beam is a convention, not geometry; the choice changes $y$ by $\bar l$ ratios
(16.67 vs 10 µm for 100 × 100 × 10 µm).

### 14.3 Tissue equivalence: $\kappa_{\rm t}$ from NIST

$$
\kappa_{\rm t}(E) = \frac{L_{\rm tissue}(E)}{L_{\rm SiC}(E)} = \frac{(S/\rho)_{\rm t}\,\rho_{\rm t}}{(S/\rho)_{\rm SiC}\,\rho_{\rm SiC}},\qquad y_{\rm tissue} = \kappa_{\rm t}(E)\,y_{\rm SiC},
$$

(`kappa_table`, `etna/core/stopping_power.py`; `compute_kappa_table(source="bragg")` and
`tissue_equivalence_correction`, `etna/core/microdosimetry.py`), with the electronic stopping
powers of NIST PSTAR/ASTAR and SiC by Bragg additivity (8.2), at the particle's **kinetic**
energy (not the energy deposited in the thin SV; audit C-2): a scalar for a mono-energetic beam
or the per-event primary energy from MC truth. There is no silent default: without a table the
caller must pass a constant. This is the **linear** (LET) ratio, below 1 because SiC is 3.21 times
denser than water; the **mass** ratio $(S/\rho)_{\rm t}/(S/\rho)_{\rm SiC}$ (~1.2-1.5, `ratio="mass"`)
converts absorbed *dose* (Bragg-Gray), not lineal energy. Earlier claims of an "inverted"
$\kappa>1$ mixed the two. **Status:** **literature data + assumption** (Bragg additivity, and
$y\propto$ LET for crossers); the equivalence of chord distributions between SiC and the tissue
site is a further assumption.

**Worked example** (`worked_examples.py microdosimetry`): $\bar l$(100 × 100 × 10 µm) = 16.6667 µm
in the code and by $4V/S$, slab 20 µm; $\kappa_{\rm t}$(proton, SiC→water) = 0.4261 / 0.3909 /
0.3784 / 0.3751 at 1 / 10 / 62 / 150 MeV; for 20 000 synthetic deposits (gamma-distributed,
mean 10 keV) on $\bar l$ = 10 µm, `lineal_moments` gives $y_F = 0.9952 \pm 0.0035$ and
$y_D = 1.2471\pm0.0054$ keV/µm, equal to the hand formulas; the shot-noise floor of 18 pA with a
1 µs shaping time (`estimate_noise_floor`, $Q_{\min} = 3\sqrt{2qI\tau}$, an intrinsic lower
bound without electronics) is 0.0072 fC, 0.351 keV, $y_{\min}$ = 0.021 keV/µm.

## 15. Monte Carlo energy deposition and straggling

### 15.1 Condensed-history sampler

`etna/core/mc_transport.py` samples, event by event, the energy a proton, alpha or C-12 deposits
in a rectangular SV (optionally inside a SiC surround), with per-event counter-based random
streams (SplitMix64 seeded from (seed, event)), so results do not depend on the thread count.
Nothing is fitted ([USAGE §11](USAGE.md#11-microdosimetric-sweep), `lab/mc_transport/notes/LOG.md`).

- **Steps.** $\Delta s = \min(s_{\rm boundary},\ f\,E/L_{\rm tot})$, $f$ = 0.05 by default; a
  particle whose CSDA range is shorter than the step deposits the rest and stops at the range.
- **Mean loss.** From the CSDA range table: $\bar\Delta_{\rm tot} = E - E(R-\Delta s)$, the
  electronic share at the mid-step energy, $\bar\Delta_e = \bar\Delta_{\rm tot}\,L_e/(L_e+L_n)$
  (NIST tables, SiC by Bragg additivity). The nuclear loss slows the ion but is not deposited as
  ionisation (no Lindhard partition).
- **Effective charge.** Protons $z=1$; alpha and C-12 the Ziegler effective charge of Geant4
  11.3 `G4ionEffectiveCharge` (C-12 with the Fermi velocity set to 1 $v_{\rm Bohr}$); C-12 stopping
  is $z_{\rm eff}^2\times$ PSTAR at the same velocity (**modelled**).
- **Delta rays above a cut $T_{\rm cut}$** (Geant4 scheme; `auto_tcut_keV`: the electron whose
  CSDA range is 5 % of the smallest SV dimension, clamped to 1-1000 keV). With
  $T_{\max} = 2m_ec^2\beta^2\gamma^2/(1+2\gamma m_e/M+(m_e/M)^2)$ and $k = 2\pi r_e^2m_ec^2N_ez^2/\beta^2$:

$$
\bar\Delta_{\rm restricted} = \bar\Delta_e - k\left[\ln\frac{T_{\max}}{T_{\rm cut}}-\beta^2\Big(1-\frac{T_{\rm cut}}{T_{\max}}\Big)\right]\Delta s,\qquad
\langle N_\delta\rangle = k\left[\Big(\frac{1}{T_{\rm cut}}-\frac{1}{T_{\max}}\Big)-\frac{\beta^2}{T_{\max}}\ln\frac{T_{\max}}{T_{\rm cut}}\right]\Delta s ,
$$

  $N_\delta$ Poisson; each $T$ sampled from $1/T^2$ on $[T_{\rm cut},T_{\max}]$ and accepted with
  probability $1-\beta^2T/T_{\max}$, emitted at the two-body angle
  $\cos^2\theta = T(T_{\max}+2m_ec^2)/(T_{\max}(T+2m_ec^2))$, and transported in a **straight line**
  over its NIST ESTAR CSDA range (`etna/core/electron_range.py`; SiC by Bragg additivity, range
  started at 1 keV with a power-law $R = T/(S(1+p))$), so delta-ray energy can leave or enter the SV.
- **Fluctuations of the restricted loss.** A line-by-line numba port of Geant4 11.3
  `G4UniversalFluctuation` / `G4UrbanFluctuation` and the Gaussian branch of `G4IonFluctuations`
  (Urban's model, Geant4 defaults) [G4]: thin layers in the Landau/Vavilov-like regime, thick ones
  Gaussian with the Bohr variance of the restricted spectrum,
  $\sigma^2 = (T_{\max}/\beta^2-T_{\rm cut}/2)\,2\pi r_e^2m_ec^2\,z^2N_e\,\Delta s$ (when
  $\bar\Delta\ge10T_{\rm cut}$ and $T_{\max}\le2T_{\rm cut}$).

**Not modelled, deliberately:** nuclear interactions and fragmentation (C-12 reaction rates of
~$10^{-4}$ per 10 µm per the Geant4 cross-section code), multiple Coulomb scattering of the primary
(straight tracks), electron scattering (straight delta rays: an upper bound on escape),
the charge-exchange factor of `G4IonFluctuations`. **Status:** mean loss **literature data**; the
fluctuation model is Geant4's (**literature model**); the sampler against the independent analytic
references of 15.2 is **verified** (`tests/test_mc_transport.py`, `tests/test_straggling.py`);
spectra of the real detector **predicted**.

### 15.2 Analytic straggling references

`etna/core/straggling.py`, independent of the sampler and used to benchmark it (free-electron
theories, no atomic binding): the Landau density with the PDG most probable value $\Delta_p$ and
width $\xi$, $\lambda = (\Delta-\Delta_p)/\xi-0.22278$ (mode at $-0.22278$, FWHM $4.018\,\xi$)
[PDG24]; the Vavilov density in the CERNLIB G116 integral form

$$
\phi(\lambda) = \frac1\pi e^{\kappa_V(1+\beta^2C)}\int_0^\infty e^{\kappa_Vf_1(u)}\cos\big(\lambda u+\kappa_Vf_2(u)\big)\,du,
$$

$f_1 = \beta^2(\ln u-\mathrm{Ci}\,u)-\cos u-u\,\mathrm{Si}\,u$,
$f_2 = u(\ln u-\mathrm{Ci}\,u)+\sin u+\beta^2\mathrm{Si}\,u$, $C$ Euler's constant,
$\lambda = (\Delta-\langle\Delta\rangle)/T_{\max}-\kappa_V(1+\beta^2-C)$ (validated by its first
three cumulants $k_n = \xi T_{\max}^{n-1}(1/(n-1)-\beta^2/n)$ to $<10^{-6}$); and the Bohr Gaussian
$\Omega_B^2 = 4\pi z^2e^4N_e\,t$. The regimes: Vavilov for $\kappa_V = \xi/T_{\max}\gtrsim0.1$,
Gaussian for $\kappa_V\gtrsim10$; for very thin absorbers ($\kappa_V\ll0.01$) real (Bichsel)
distributions are wider than Landau's, as the PDG states.

**Worked example** (`worked_examples.py microdosimetry`): $N_e$(SiC) = $9.6423\times10^{23}$
cm$^{-3}$ and the Bohr width of a proton crossing 10 µm of SiC,
$\Omega_B = \sqrt{4\pi e^4N_et}$ = 15.85 keV; the Bohr range straggling of a 62 MeV proton in SiC
is 147.5 µm (`worked_examples.py generation`).

## 16. Numerics

### 16.1 Nonlinear solvers and tolerances

- **devsim** (`etna/core/drift_diffusion.py`, `etna/core/solver_tol.py`): Newton on the coupled
  system with devsim's `log_damp` / `positive` updates. A DC solve is accepted at relative update
  $10^{-8}$ and absolute $10^{10}$ (strict tier); fallback tiers are relaxed ($10^{12}$ /
  $10^{-8}$) and aggressive ($10^{14}$ / $10^{-6}$, 200 iterations). A solve that hits its
  iteration cap with its final update already inside the aggressive tolerance is accepted as an
  aggressive-tier point (`solve_tier`: the limit cycle of zero-flux surfaces next to lateral p$^+$
  implant edges, RelError ~$5\times10^{-8}$). Every fallback is logged and counted
  (`ramp_relaxed_retries`, `ramp_aggressive_retries`). On 1D devices $10^{-8}$ gives outputs
  bit-identical to $10^{-10}$ (Newton is quadratic); the FLASH continuation keeps $10^{-10}$.
- **Batched Poisson**: per-device Newton with a Thomas solve, update clamped at 20 $V_t$, stop when
  the max update is below the tolerance; a secant predictor between requested biases (no fixed
  0.5 V lattice).
- **Gummel** (batched and 2D engines): electrons with $p$ frozen in $U$, then holes, then the
  quasi-Fermi Poisson; in 2D with exact residuals and the three dark-state band LUs as chord
  corrections, so the fixed point is the exact discrete solution. With avalanche the inner
  $(n,p)$ problem is a full Newton (2×2-block tridiagonal) at frozen $\psi$ (7.2).
- **Track engine**: fully coupled Newton on interleaved $(\delta\psi,n,p)$ per node, a modified
  (chord) Newton that reuses one SuperLU factorisation; the excess form keeps the edge potential
  differences at full precision in the $10^{19}$ layers (`etna/core/track/kernels.py`).
- **Charged traps**: Newton on $\psi$ with the analytic trap-charge derivative, clamp 2 $V_t$, up to
  2000 iterations (`TRAP_CLAMP`, `TRAP_MAX_ITER`).

### 16.2 Bias stepping and seeded bias

`ramp_bias` / `_ramp_bias_adaptive`: steps between 0.5 and 2 V, halved on a strict failure,
doubled after success, the Newton start for $\psi$ linearly extrapolated from the last two
converged states. `seeded_bias` replaces the ramp by an **engine seed** (the batched / 2D dark
state at the target bias, carriers floored) followed by **one strict devsim solve**; it is accepted
only if that solve converges strictly, otherwise the solution is restored and the ramp runs, so the
fallback is the old behaviour bit for bit (`tests/test_seeded_bias.py`). The seed is skipped for
forward bias, laterally non-uniform or cylindrical 2D devices, generation present, and physics the
engine lacks (Auger, TAT, opt-in ionization) unless a call site allows it.

### 16.3 Time integration

- devsim transients: BDF1 with the step rules of 9.3 and no local error control.
- Track engine (`etna/core/track/solver.py`): TR-BDF2, $\gamma = 2-\sqrt2$, a trapezoidal stage to
  $t+\gamma h$ and a BDF2 stage to $t+h$ (L-stable, one-step), with the embedded error estimate of
  Bank et al. 1985 as given in the module docstring (the paper was not reopened by this project):

$$
\mathrm{LTE}\approx2Ch\left[\frac{f_n}{\gamma}-\frac{f_\gamma}{\gamma(1-\gamma)}+\frac{f_{n+1}}{1-\gamma}\right],\qquad C = \frac{-3\gamma^2+4\gamma-2}{12(2-\gamma)},
$$

  on the carrier rates; the step is accepted when the volume-integrated $|\mathrm{LTE}|$ of both
  carriers is below `rtol` (default $10^{-4}$) times the generated pairs, and resized by
  $0.9E^{-1/3}$. The step is an exact quadrature of the discrete rates with weights
  $w_0 = w_1 = h/(2(2-\gamma))$, $w_2 = h(1-\gamma)/(2-\gamma)$, which closes the pair ledger.
- Trap kinetics: backward Euler in $f$ with step doubling (12.5).

### 16.4 Meshes and adaptivity

- 1D: devsim's graded line mesh (1.1); the engines read its nodes. `run_breakdown` re-solves on
  the bisected mesh for `mesh_error_V` (7.2).
- 2D builders: tensor meshes; `etna/core/contact_mesh.py` grades the contact-on-epi cathode
  (`depth_mesh="graded"`, default since 2026-10-21: 0.1 nm cathode cell growing by 1.25$\times$ to
  50 nm, 1 nm at the junction; `"legacy"` = 100 nm / 10 nm).
- fv2d adaptive refinement (`etna/core/fv2d/refine.py`, opt-in): tensor-line bisection on high
  potential steps, high curvature, and the contact cells; goal-oriented stop on the CCE scan with
  the geometric extrapolation $d_L\,r/(1-r)$, $r = d_L/d_{L-1}$ clipped to [0.1, 0.9]. **The
  estimate is not a bound** (on the mesa centre at 10 V it stopped at $7.8\times10^{-6}$ with a true
  error of $1.2\times10^{-4}$); for reference values use uniform bisection with Richardson
  extrapolation once increments are asymptotic.
- Track engine: first radial cell 0.5-2 nm, geometric growth to $R$ = 20-30 µm; `MESH_FAST` PHD
  values carry 5-8 % mesh error (USAGE §19).

### 16.5 The contact-cell finding

An ohmic contact on the n$^-$ epi is a perfect sink for holes. Holes generated in its boundary
layer (the undepleted epi at partial depletion, or a drift-diffusion layer of thickness
$\sim k_BT/(qE)$, 13-26 nm at 1-2$\times10^4$ V/cm, once the depletion reaches the contact) are
lost. The box scheme lumps the contact cell's generation onto the contact node, where the contact
equation replaces continuity (4.2), so the CCE converges at **first order in the contact-cell size**
and a coarse cell **overestimates** it: +6.2e-3 at 10 V and +7.9e-4 at 50 V on the 100 nm cell of
the alternative-structure builders (mesa centre, against the converged 0.9630651 / 0.9990906);
with the graded mesh the error is +5.4e-5 / $-2.9\times10^{-6}$ (`etna/core/contact_mesh.py`,
`lab/mesh_contact/notes/LOG.md`). Field and curvature indicators alone never refine that cell and
stop with a false "converged". The default n$^+$ stack does not have this contact: its cathode sits
on the $10^{19}$ layer, generation is confined to the epi and the high-low junction reflects holes;
there the 1D engine moves by $<2\times10^{-10}$ when its contact cells are bisected.

A related artefact: with the cathode on the epi (legacy stacks), the field in the last cell grows
with bias after reach-through and devsim stops converging when it reaches ~$4.7\times10^4$ V/cm
(1D near 60-66 V for the legacy device; 2D 54.5-56.5 V for $N_b$ 5-8.5$\times10^{13}$ on the 10 µm epi, `etna/core/engine2d/scan.py`). It is a
contact artefact, **not** punch-through; ETNA reports the truncation (`truncated`, NaN rows,
`Engine2DOutsideValidity`) instead of extrapolating. The n$^+$ layer removes it.

### 16.6 Code-to-code verification

Inside documented parameter boxes (USAGE §7) the engines agree with devsim to:
$|\Delta C|/C\lesssim3\times10^{-8}$ and $|\Delta\mathrm{CCE}|\lesssim10^{-9}$ in 1D to $-100$ V,
2D potential ~$10^{-13}$ V, FLASH CCE to $\le7.5\times10^{-12}$ with Auger, and the dark current within its own parity tests
(`tests/test_batched_engine.py`, `tests/test_engine2d*.py`, `tests/test_dark_current_engine.py`,
`tests/test_flash_engine.py`, `tests/test_fv2d.py`). These numbers show that the equations are
solved consistently, not that they describe the device.

**Worked example** (`worked_examples.py numerics`): the small-signal $C$ of the default device from
the batched engine and from devsim's ±1 mV difference agree to $9.0\times10^{-9}$ / $1.6\times10^{-9}$ /
$6.9\times10^{-11}$ relative at 0 / $-5$ / $-10$ V ($C(0) = 2.07608\times10^{-10}$ F on 4 mm$^2$).

## 17. Uncertainty

ETNA has **no end-to-end uncertainty propagation**: no output carries an error bar that combines
parameter, model and numerical uncertainty. What exists:

| kind | what | where |
| --- | --- | --- |
| numerical: mesh | `mesh_error_V` of $V_{BD}$ (bisected-mesh re-solve), `quadrature_error_V` | `etna/api/models.py`, `etna/core/avalanche.py` |
| numerical: mesh | fv2d refinement estimate (not a bound, 16.4); contact-cell convergence studies | `etna/core/fv2d/refine.py`, `lab/mesh_contact` |
| numerical: time | TR-BDF2 LTE control, pair ledger residual, energy-theorem residual | `etna/core/track/solver.py` |
| numerical: solver | tolerance tiers counted per sweep; truncation flags | `etna/core/drift_diffusion.py` |
| statistical | standard errors of $y_F$, $y_D$ (delta method), adaptive sampling to a target | `etna/core/mc_response.py` |
| parameter: fit | Gauss-Newton covariance (log space), reduced-$\chi^2$ scaling, correlations, singular values and identifiability flags of trap-rate scales fitted to C-V | `etna/core/trap_calibration.py` |
| parameter: fit | the ± of the calibrated direct removal (data spread over depletion widths) | `etna/core/trap_presets.py` |
| parameter: bracket | `cce_uncertainty_envelope`: min/max CCE over the 8 corners of $\times0.5$/$\times2$ on the three introduction rates | `etna/core/radiation_damage.py` |
| parameter: alternatives | registries of literature sets (mobility, $v_{\rm sat}$, $\alpha$) whose spread is a model-uncertainty proxy, e.g. the $V_{BD}$ spread across impact-ionisation sets | `etna/core/physics_models.py`, USAGE §18 |
| model: stated | provenance notes with ranges: electron mobility ±8 %, holes $-35/+0$ % along c, Auger $10^{-31}$-$2\times10^{-30}$, EH6/7 $\sigma_n$ $6\times10^{-14}$-$2\times10^{-13}$, $r_n$ 1.02-1.20, Bragg additivity a few % (unverified), Hurkx closed form ~30 % in shape | the cited modules |
| data | plausibility checks of measurement files (unit mix-ups, impossible values) | `etna/measurements/plausibility.py` |

The dominant uncertainties of the physical predictions are not numerical: the doping profile is
inferred from one C-V curve without SIMS, the lifetimes of this epi are unmeasured (a tilted-alpha
fit wants ~23 ns, not the library 1 ns / 600 ns), the leakage shape is a heuristic, the damage
removal mechanism is unidentified and calibrated on other devices, and the SiC neutron NIEL
normalisation is computed but unchecked. The numerical errors of the default outputs (parity
$\le10^{-8}$, mesh errors $\le10^{-4}$ in CCE on the graded meshes) are orders of magnitude below
these.

## References

Only sources that this project's agents recorded as opened (in the GitHub issues #3-#9, `lab/*/notes`,
`.planning/` or the data READMEs) are listed as references; "via X" means the numbers were read in
the opened source X, not in the primary. Citations that appear in the code but were never opened are
listed separately and must be treated as unverified.

**Opened**

- [AYA] T. Ayalew, *SiC semiconductor devices technology, modeling and simulation*, PhD thesis, TU Wien
  (2004), https://www.iue.tuwien.ac.at/phd/ayalew/ (nodes 61-75: mobility Tables 3.4-3.6, eqs.
  3.72-3.81, Sec. 3.5.3, Sec. 3.6).
- [BUR24r] J. Burin et al., "TCAD parameters for 4H-SiC: a review", arXiv:2410.06798 (v8, 2025). Via it:
  Schaffer 1994, Lades 2000, Nilsson 2000, Hatakeyama 2003 (mobility, $v_{\rm sat}$, anisotropy);
  Niwa 2015, Zhao 2019, Loh 2008, Stefanakis 2021 (impact ionisation); Kagamihara 2004, Achatz 2008
  (ionization); Galeckas et al., APL 71, 3269 (1997) [GAL97] and Tanaka 2023 (Auger); Lades
  (carrier-carrier scattering).
- [BUR24t] J. Burin et al., arXiv:2407.16710 (2024), Table I (trap levels, cross-sections,
  introduction rates).
- [BV23] M. Bruzzi, E. Verroi, "Epitaxial SiC dosimeters and flux monitoring detectors for proton
  therapy beams", Materials 16 (2023) 3643, doi:10.3390/ma16103643 (read via PMC10223321).
- [CAP18] I. Capan et al., J. Appl. Phys. 124, 245701 (2018), arXiv:1812.06462, Table I.
- [FRE02] E. Fretwurst, G. Lindström, I. Pintilie, J. Stahl, arXiv:physics/0211118 (the 95 MeV mb
  normalisation).
- [G4] Geant4 11.3 source: `G4UniversalFluctuation`, `G4UrbanFluctuation`, `G4IonFluctuations`,
  `G4ionEffectiveCharge` (L. Urban's fluctuation model, Geant4 Physics Reference Manual).
- [GSP23] A. Gsponer et al., arXiv:2310.02047 (TRIGA neutron irradiation of 4H-SiC, neq labelling).
- [HAT04] T. Hatakeyama, T. Watanabe, T. Shinohe, K. Kojima, K. Arai, N. Sano, Appl. Phys. Lett. 85,
  1380 (2004), doi:10.1063/1.1784520, Table I; and T. Hatakeyama, SISPAD 2005, pp. 171-174, Eq. 1.
- [IOFFE] Ioffe Institute NSM archive, SiC, https://www.ioffe.ru/SVA/NSM/Semicond/SiC/ (opened for the
  thermal values $c_p$ and diffusivity only; its other values in `sic_material.py` are cited in the
  code but were not re-checked).
- [ISH25] R. Ishikawa, *Low-field and high-field anisotropic carrier transport in SiC*, PhD thesis,
  Kyoto University (2025), doi:10.14989/doctor.k26017, http://hdl.handle.net/2433/296110 (Eqs. 2.16,
  2.17, 4.5, 4.6, Table 3.1, Sec. 4.3.1); cross-check R. Ishikawa, M. Kaneko, T. Kimoto, Solid State
  Phenom. 360, 205 (2024), doi:10.4028/p-MAe9pZ.
- [KOR14] Korcyl, arXiv:1410.5250, eqs. 1.20 and 2.3-2.12 (via it: Waligórski, Hamm & Katz 1986;
  Zhang et al. 1985).
- [NIST] M. J. Berger, J. S. Coursey, M. A. Zucker, J. Chang, ESTAR, PSTAR and ASTAR, NIST,
  https://physics.nist.gov/PhysRefData/Star/Text/ (tables downloaded 2026-10-20/21; URLs and SHA-256 in
  `data/nist_star/README.md`, `data/nist_estar/README.md`).
- [P25a] G. Petringa et al., Physica Medica 140 (2025) 105675 (microdosimeter geometry, pp. 2-4; its
  Fig. 4b is the tilted-alpha CCE vs bias, not the C-V figure below).
- [P25b] G. Petringa et al., JINST 20 (2025) C08019.
- [PDG24] Particle Data Group, Review of Particle Physics 2024, "Passage of particles through matter"
  (via it: Vavilov 1957, Bichsel 1988).
- [PHOT] The "Photons" paper of the INFN-LNS group, held in the repository as
  `papers/SiC_Photons_MedicalPhysics (3).pdf`; its Fig. 4b (C-V depletion width, PDF page 10) is the
  calibration curve of the default doping. **No bibliographic record (journal, volume, DOI) is stored
  in the repository**; do not confuse it with [P25a].
- [QIA17] X. Qian, P. Jiang, R. Yang, arXiv:1712.00830, Table 1 (thermal conductivity of 4H-SiC).
- [RAD25] Radmanovac et al., arXiv:2510.11304, Table 2 (252.7 MeV proton C-V slopes).
- [SRNIEL] SR-NIEL web calculator, Protons & Ions NIEL Calculator, SR-NIEL-7 ver. 11.1 (12 May 2026),
  https://www.sr-niel.org/ (fetched 2026-10-20; raw pages in `data/niel/`).
- NJOY2016.78 HEATR on ENDF/B-VIII.0, computed in this repository (`scripts/build_neutron_niel.py`,
  `data/niel/README.md`); Griffin et al. SAND92-0094 damage function (RD48 compilation) stored
  verbatim for the Si check.
- DLTS checkpoints of EH6/7 as opened in `lab/radiation_traps/notes/MECHANISMS.md`: Hazdra et al.,
  SSP 242 (2016) 421; Danno & Kimoto, JAP 100 (2006) 113728; Alfieri & Kimoto, APL 102 (2013) 152108.

**Cited in the code but not opened by this project (unverified)**

Hurkx, Klaassen & Knuvers, IEEE TED 39, 331 (1992) (the $\Gamma$ closed form is derived in the
`hurkx_gamma` docstring, not checked against the paper); Scharfetter & Gummel (1969); McIntyre (1966);
Bank et al., IEEE TED 32, 1992 (1985) (TR-BDF2 error estimate); Hecht, Z. Phys. 77, 235 (1932); ICRU
Report 36 (1983); Hazdra et al., NIM B 327 (2014) (origin of the Table I rates; title only); Pensl &
Choyke (1993) and Kimoto & Cooper (legacy ionization); "PRB 86, 075205" (annealing barrier, no authors
recorded); the "Lebedev" W-value of 7.7-7.8 eV; arXiv:2408.12410 (proton hardness systematic, quoted
without an opened marker); the NIM A literature fit behind the 7 µm hole diffusion length (PII
S0168900205006443, authors not recorded).


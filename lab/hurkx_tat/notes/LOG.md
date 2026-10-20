# Hurkx TAT enhancement (audit-v7 B-08), 2026-10-20

Commands: `uv run python lab/hurkx_tat/bench/field_gamma.py`,
`uv run python lab/hurkx_tat/bench/nt_refit.py` (run on the pre-refit N_t; nt_refit.py
rescales N_t linearly to I(300 K, -60 V, 4 mm2) = 18 pA, devsim, 31-point sweep).

## field_gamma.py (default device, devsim, old model, Gamma closed form evaluated post hoc)

    F_Gamma 1.304e5 / 3.688e5 / 6.775e5 V/cm at 150 / 300 / 450 K (m_t = 0.25)
    -10 V  Emax 1.711e5 V/cm  Gamma_max 3.53  mean Gamma (E > E_ref) 1.18
    -60 V  Emax 1.711e5       Gamma_max 3.53  mean 1.28
    -100 V Emax 1.761e5       Gamma_max 3.68  mean 1.97
    I_SRH ~ -4e-16 A (~1e-5 of I_total): the SRH share of the enhancement is negligible

## nt_refit.py (Hurkx model in place)

    default:          N_t 2.849024e12  I60 37.8308 pA -> N_t 1.355575e12  I60 18.000000 pA
    tabulated:        N_t 2.849219e12  I60 37.8182 pA -> N_t 1.356120e12  I60 18.000000 pA
    nplus_calibrated: N_t 2.805336e12  I60 37.5723 pA -> N_t 1.343970e12  I60 18.000000 pA
    legacy:           N_t 2.800000e12  I60 37.3735 pA -> (1.348550e12, NOT applied: the
                      legacy preset keeps its historical N_t)
    microdosimeter:   area 0.25 cm2, uses the default N_t (237.9 pA -> not a calibrated preset)

## Closed form vs Hurkx integral a int_0^1 exp(a u - K u^1.5) du (Delta E = E_t = 0.67 eV)

Gamma (integral / closed form):

    150 K: 5e4 V/cm 2.120 / 2.728; 1e5 7.417 / 8.480; 1.7e5 42.04 / 43.82
    300 K: 5e4 0.653 / 0.848; 1e5 1.355 / 1.792; 1.7e5 2.788 / 3.501; 3e5 8.568 / 9.681
    450 K: 5e4 0.381 / 0.456; 1e5 0.708 / 0.926; 1.7e5 1.238 / 1.641

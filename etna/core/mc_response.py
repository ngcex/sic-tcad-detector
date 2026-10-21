"""Event-by-event detector response of sampled tracks: lineal energy with adaptive
statistics, and pulse heights through the TCAD CCE(LET) table.

Pipeline (each step opt-in, nothing here changes the existing ``process_mc_ensemble``
or ``run_microdosimetry`` behaviour):

    sample_track_events  ->  epsilon_k, chord_k  (etna.core.mc_transport)
    pulse_heights        ->  N_pairs (W, optional Fano) -> x CCE(LET_k) -> + ENC noise
                             -> threshold -> measured energy per event
    lineal_moments / adaptive_lineal_energy -> y_F, y_D with standard errors

LET for the CCE lookup is per event: ``let_mode="chord"`` uses epsilon_k / chord_k
(the event's own path; for an oblique or isotropic track this is the deposited energy
per unit path, not epsilon / thickness), ``"thickness"`` the old epsilon / t convention
of ``process_mc_ensemble``.
"""

from __future__ import annotations

import math

import numpy as np

from etna.core.sic_material import SiC4H_Parameters

E_PAIR_SIC_EV = SiC4H_Parameters().E_pair_eV


def lineal_moments(edep_keV, l_bar_um):
    """y_F, y_D (keV/um) of the events with epsilon > 0 and their standard errors.

    Event-level estimators (no binning): y_F = <eps>/l, y_D = <eps^2>/(<eps> l). The
    standard error of the ratio y_D uses the delta method with the sample covariance
    of (eps, eps^2).
    """
    e = np.asarray(edep_keV, dtype=float)
    e = e[e > 0]
    n = e.size
    if n < 2:
        return {
            "y_F": np.nan,
            "y_D": np.nan,
            "se_y_F": np.nan,
            "se_y_D": np.nan,
            "n": n,
        }
    y = e / l_bar_um
    m1 = y.mean()
    y2 = y * y
    m2 = y2.mean()
    v1 = y.var(ddof=1)
    v2 = y2.var(ddof=1)
    c12 = np.cov(y, y2, ddof=1)[0, 1]
    yD = m2 / m1
    # d(yD) = dm2/m1 - m2 dm1/m1^2
    var_yD = (v2 / m1**2 + (m2**2 / m1**4) * v1 - 2.0 * (m2 / m1**3) * c12) / n
    return {
        "y_F": float(m1),
        "y_D": float(yD),
        "se_y_F": float(math.sqrt(v1 / n)),
        "se_y_D": float(math.sqrt(max(var_yD, 0.0))),
        "n": int(n),
    }


def default_l_bar_um(sv_um, irradiation):
    """Mean chord of the SV for the sampler's irradiation (``sample_track_events``).

    Isotropic (mu-randomness): Cauchy 4V/S. Parallel beam tilted by theta in the x-z
    plane, spread uniformly over the projected area: V / A_proj with
    A_proj = wx wy cos(theta) + wy t sin(theta) (mean chord of a convex body for a
    uniform parallel beam). This is NOT the infinite-slab path t / cos(theta), which
    overestimates it once the side face is lit (60 deg, 100x100x10 um: 17.05 vs
    20 um; issue #18). Normal beam: t.
    """
    wx, wy, t = (float(v) for v in sv_um)
    V = wx * wy * t
    if irradiation == "isotropic":
        S = 2.0 * (wx * wy + wx * t + wy * t)
        return 4.0 * V / S
    if irradiation == "normal":
        return t
    th = math.radians(float(irradiation))
    return V / (wx * wy * math.cos(th) + wy * t * math.sin(th))


def adaptive_lineal_energy(
    particle,
    energy_MeV,
    rel_tol=0.01,
    batch=50000,
    max_events=5_000_000,
    l_bar_um=None,
    seed=0,
    keep_events=True,
    **sample_kw,
):
    """Sample batches until the relative standard errors of y_F AND y_D are < rel_tol.

    ``sample_kw`` go to :func:`etna.core.mc_transport.sample_track_events` (sv_um,
    material, irradiation, tcut_keV, ...). Batches use consecutive event indices of
    one random stream (``event_offset``), so the result for a given (seed, number of
    events) does not depend on the batch size. Returns y_F, y_D, their standard errors,
    the number of events, ``converged`` and the per-batch ``history``.
    """
    from etna.core.mc_transport import sample_track_events

    if not (0 < rel_tol < 1):
        raise ValueError("rel_tol must be in (0, 1)")
    sv = sample_kw.get("sv_um", (100.0, 100.0, 10.0))
    irr = sample_kw.get("irradiation", "normal")
    lb = default_l_bar_um(sv, irr) if l_bar_um is None else float(l_bar_um)
    chunks, chords, hist = [], [], []
    n = 0
    mom = None
    converged = False
    meta = None
    while n < max_events:
        k = int(min(batch, max_events - n))
        ev = sample_track_events(
            particle, energy_MeV, n_events=k, seed=seed, event_offset=n, **sample_kw
        )
        meta = ev.meta
        chunks.append(ev.edep_keV)
        chords.append(ev.chord_um)
        n += k
        mom = lineal_moments(np.concatenate(chunks), lb)
        rF = mom["se_y_F"] / mom["y_F"] if mom["y_F"] > 0 else np.inf
        rD = mom["se_y_D"] / mom["y_D"] if mom["y_D"] > 0 else np.inf
        hist.append({"n_events": n, "rel_se_y_F": rF, "rel_se_y_D": rD, **mom})
        if rF < rel_tol and rD < rel_tol and len(hist) >= 2:
            converged = True
            break
    out = {
        **mom,
        "l_bar_um": lb,
        "n_events": n,
        "converged": converged,
        "rel_tol": rel_tol,
        "history": hist,
        "meta": meta,
    }
    if keep_events:
        out["edep_keV"] = np.concatenate(chunks)
        out["chord_um"] = np.concatenate(chords)
    return out


def pulse_heights(
    edep_keV,
    chord_um=None,
    cce_interp=None,
    let_mode="chord",
    sv_thickness_um=None,
    W_eV=E_PAIR_SIC_EV,
    fano=None,
    enc_electrons=0.0,
    threshold_keV=0.0,
    seed=0,
):
    """Event-by-event pulse heights from deposited energies.

    N = eps / W pairs (Gaussian with variance F N when ``fano`` is given; None = no
    pair statistics), Q = CCE(LET) N with ``cce_interp`` from
    ``single_particle.load_cce_let_table`` (None: CCE = 1), plus Gaussian electronic
    noise of ``enc_electrons`` (rms e-), then the threshold on the measured energy
    Q W. Returns arrays per INPUT event (eps <= 0 gives LET 0 and Q from noise only)
    and ``triggered``.

    The Allpix Squared default for 4H-SiC is W = 7.6 eV, F = 0.1 (Allpix user manual,
    "Sensor Material Properties"); ETNA's W is 7.8 eV (SiC4H_Parameters). F is left
    to the caller (no default) because no value has been checked against SiC data here.
    """
    e = np.asarray(edep_keV, dtype=float)
    if np.any(~np.isfinite(e)) or np.any(e < 0):
        raise ValueError("edep_keV must be finite and >= 0")
    if let_mode not in ("chord", "thickness"):
        raise ValueError("let_mode must be 'chord' or 'thickness'")
    if W_eV <= 0 or enc_electrons < 0 or threshold_keV < 0:
        raise ValueError("W_eV > 0, enc_electrons >= 0, threshold_keV >= 0 required")
    if fano is not None and not (0 <= fano <= 1):
        raise ValueError("fano must be in [0, 1]")
    if let_mode == "chord":
        if chord_um is None:
            raise ValueError("let_mode='chord' needs chord_um")
        c = np.asarray(chord_um, dtype=float)
        if sv_thickness_um is not None:
            c = np.where(c > 0, c, float(sv_thickness_um))
        let = np.where(c > 0, e / np.where(c > 0, c, 1.0), 0.0)
    else:
        if sv_thickness_um is None:
            raise ValueError("let_mode='thickness' needs sv_thickness_um")
        let = e / float(sv_thickness_um)
    rng = np.random.default_rng(seed)
    n_mean = e * 1e3 / W_eV
    if fano:
        n_pairs = np.maximum(
            n_mean + rng.standard_normal(e.size) * np.sqrt(fano * n_mean), 0.0
        )
    else:
        n_pairs = n_mean
    if cce_interp is None:
        cce = np.ones_like(e)
        n_clamped = 0
    else:
        pos = let > 0
        cce = np.ones_like(e)
        if np.any(pos):
            vec = getattr(cce_interp, "vectorized", False)
            vals = (
                np.asarray(cce_interp(let[pos]), dtype=float)
                if vec
                else np.array([cce_interp(v) for v in let[pos]])
            )
            cce[pos] = vals
        lo = getattr(cce_interp, "let_min", None)
        hi = getattr(cce_interp, "let_max", None)
        n_clamped = (
            int(np.sum(pos & ((let < lo) | (let > hi)))) if lo is not None else 0
        )
    q = cce * n_pairs
    if enc_electrons > 0:
        q = q + rng.standard_normal(e.size) * enc_electrons
    measured = q * W_eV * 1e-3
    trig = measured > threshold_keV
    return {
        "LET_keV_um": let,
        "n_pairs": n_pairs,
        "CCE": cce,
        "charge_e": q,
        "charge_fC": q * 1.602176634e-4,
        "measured_keV": measured,
        "triggered": trig,
        "n_triggered": int(trig.sum()),
        "n_clamped": n_clamped,
        "W_eV": float(W_eV),
        "fano": fano,
        "enc_electrons": float(enc_electrons),
        "threshold_keV": float(threshold_keV),
    }

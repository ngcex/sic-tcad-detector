"""Figures of the measurement comparisons in the etna.viz house style.

Each function returns a matplotlib ``Figure`` built under ``etna.viz.style.RC`` with a
provenance footer (data file + hash, model solver + config hash, date) and the
measurement-comparison label; export with ``etna.viz.export_figure(fig, "png")``.
Nothing is simulated here: every point comes from the result passed in.
"""

from __future__ import annotations

import matplotlib as mpl
import numpy as np
from matplotlib.figure import Figure
from matplotlib.layout_engine import ConstrainedLayoutEngine

from etna.viz import style as S
from etna.measurements._common import footer, footer_line


def _fig(size, rows=2, ratios=(2.2, 1.0), sharex=True):
    fig = Figure(figsize=S.size_inches(size), layout=ConstrainedLayoutEngine())
    if rows == 1:
        return fig, [fig.add_subplot()]
    gs = fig.add_gridspec(rows, 1, height_ratios=ratios)
    a0 = fig.add_subplot(gs[0])
    a1 = fig.add_subplot(gs[1], sharex=a0 if sharex else None)
    return fig, [a0, a1]


def cv_comparison_figure(cmp, size="single") -> Figure:
    """W = eps*A/C vs reverse bias (data and model) and the W residual."""
    with mpl.rc_context(S.RC):
        fig, (ax, ar) = _fig(size)
        v = -cmp.V
        if cmp.sigma_W_um is not None:
            ax.errorbar(
                v,
                cmp.W_data_um,
                yerr=cmp.sigma_W_um,
                fmt="o",
                ms=2.2,
                lw=0.6,
                color=S.DATA_COLOR,
                capsize=1.0,
                label="measured",
            )
        else:
            ax.plot(v, cmp.W_data_um, "o", ms=2.2, color=S.DATA_COLOR, label="measured")
        ax.plot(v, cmp.W_model_um, "-", color=S.PALETTE[0], label="model (run_cv)")
        ax.set_ylabel(r"$\varepsilon A/C$ ($\mu$m)")
        ax.set_title(
            f"C-V vs model: rms W {cmp.rms_W_um:.3g} um, rms C {100 * cmp.rms_C_rel:.3g} %"
        )
        ax.legend(loc="lower right")
        ar.axhline(0.0, color="#666666", lw=0.6)
        ar.plot(v, cmp.resid_W_um, "o-", ms=2.0, lw=0.6, color=S.PALETTE[1])
        ar.set_xlabel(S.UNITS["reverse_bias"][1])
        ar.set_ylabel(r"data $-$ model ($\mu$m)")
        footer(
            fig, cmp.notes + [footer_line(cmp.data_provenance, cmp.model_provenance)]
        )
    return fig


def inversion_figure(inv, dep=None, size="single") -> Figure:
    """Inverted donor table vs depth (probed range shaded) and the W fit residual."""
    with mpl.rc_context(S.RC):
        fig, (ax, ar) = _fig(size, ratios=(1.6, 1.0), sharex=False)
        lo, hi = inv.probed_depth_um
        ax.axvspan(lo, hi, color="#e8f1fa", lw=0, label="probed by the C-V")
        N = 10.0**inv.log10_N
        err = N * np.log(10.0) * np.nan_to_num(inv.sigma_log10_N)
        ax.errorbar(
            inv.knots_um,
            N,
            yerr=err,
            fmt="s-",
            ms=2.6,
            lw=0.9,
            color=S.PALETTE[0],
            capsize=1.0,
            label="inverted donor_table",
        )
        if dep is not None:
            ok = dep.valid
            ax.plot(
                dep.depth_um[ok],
                dep.N_cm3[ok],
                ".",
                ms=2.0,
                color=S.PALETTE[2],
                label="depletion approximation",
            )
        ax.set_yscale("log")
        ax.set_xlabel(r"Depth below the junction ($\mu$m)")
        ax.set_ylabel(r"$N_D$ (cm$^{-3}$)")
        ax.legend(loc="upper center", fontsize=S.BASE_FONT_PT - 2)
        ar.axhline(0.0, color="#666666", lw=0.6)
        ar.plot(-inv.V, inv.resid_W_um, "o-", ms=2.0, lw=0.6, color=S.PALETTE[1])
        ar.set_xlabel(S.UNITS["reverse_bias"][1])
        ar.set_ylabel(r"W residual ($\mu$m)")
        ar.set_title(
            f"fit rms {inv.rms_W_um:.3g} um, chi2 {inv.chi2:.3g} / {inv.dof} dof, lam {inv.lam:g}",
            fontsize=S.BASE_FONT_PT - 1,
        )
        footer(
            fig, inv.notes + [footer_line(inv.data_provenance, inv.model_provenance)]
        )
    return fig


def profile_figure(profile, table, config=None, size="single") -> Figure:
    """SIMS/SRP N (and Al) vs depth below the junction with the resampled knots and,
    with ``config``, the model's donor profile."""
    with mpl.rc_context(S.RC):
        fig, (ax,) = _fig(size, rows=1)
        z = profile.depth_um - table.junction_depth_um
        if profile.N is not None:
            m = profile.N > 0
            ax.plot(
                z[m],
                profile.N[m],
                "-",
                lw=0.7,
                color=S.DATA_COLOR,
                label=f"{profile.technique} N",
            )
        if profile.Al is not None:
            m = profile.Al > 0
            ax.plot(
                z[m],
                profile.Al[m],
                "-",
                lw=0.7,
                color=S.PALETTE[3],
                label=f"{profile.technique} Al",
            )
        err = table.N_D * np.nan_to_num(table.rel_uncertainty)
        ax.errorbar(
            table.knots_um,
            table.N_D,
            yerr=err,
            fmt="s",
            ms=2.6,
            color=S.PALETTE[0],
            capsize=1.0,
            label="donor_table knots",
        )
        if table.censored.any():
            ax.plot(
                table.knots_um[table.censored],
                table.N_D[table.censored],
                "v",
                ms=4,
                color=S.PALETTE[1],
                label="at detection limit",
            )
        if config is not None:
            from etna.measurements.sims import model_donor_profile

            zz = np.linspace(0.0, table.epi_thickness_um, 400)
            ax.plot(
                zz,
                model_donor_profile(config, zz),
                "--",
                lw=0.9,
                color=S.PALETTE[2],
                label="model N_D (config)",
            )
        ax.axvline(0.0, color="#666666", lw=0.6, ls=":")
        ax.axvline(table.epi_thickness_um, color="#666666", lw=0.6, ls=":")
        ax.set_yscale("log")
        ax.set_xlabel(r"Depth below the junction ($\mu$m)")
        ax.set_ylabel(r"Concentration (cm$^{-3}$)")
        ax.legend(loc="best", fontsize=S.BASE_FONT_PT - 2)
        footer(fig, table.notes + [footer_line(table.data_provenance)])
    return fig


def arrhenius_figure(analysis, size="single") -> Figure:
    """ln(|I|/T^p) vs 1/kT per bias (data, fit, model) and the ln residuals."""
    from etna.measurements._common import K_B_EV

    with mpl.rc_context(S.RC):
        fig, (ax, ar) = _fig(size)
        p = analysis.T_power
        model = {c.bias_V: c.model_fit for c in analysis.comparisons}
        for i, f in enumerate(analysis.fits):
            col = S.PALETTE[i % (len(S.PALETTE) - 1)]
            x = 1.0 / (K_B_EV * f.T_K)
            lab = "" if f.bias_V is None else f" {f.bias_V:g} V"
            ax.plot(
                x,
                f.I_abs_A / f.T_K**p,
                "o",
                ms=2.6,
                color=col,
                label=f"data{lab}: E_a {f.E_a_eV:.3f} +/- {f.sigma_E_a_eV:.3f} eV",
            )
            ax.plot(x, f.I_fit() / f.T_K**p, "-", lw=0.8, color=col)
            mf = model.get(f.bias_V)
            if mf is not None:
                ax.plot(
                    x,
                    mf.I_abs_A / f.T_K**p,
                    "--",
                    lw=0.8,
                    color=col,
                    label=f"model{lab}: E_a {mf.E_a_eV:.3f} eV",
                )
            ar.plot(x, f.residual_ln, "o-", ms=2.0, lw=0.6, color=col)
        ax.set_yscale("log")
        ax.set_ylabel(r"$|I|/T^{%g}$ (A K$^{-%g}$)" % (p, p) if p else "|I| (A)")
        ax.legend(loc="best", fontsize=S.BASE_FONT_PT - 2)
        ar.axhline(0.0, color="#666666", lw=0.6)
        ar.set_xlabel(r"$1/k_BT$ (eV$^{-1}$)")
        ar.set_ylabel("ln residual")
        footer(
            fig,
            analysis.notes
            + [
                footer_line(analysis.data_provenance, analysis.model_provenance or None)
            ],
        )
    return fig

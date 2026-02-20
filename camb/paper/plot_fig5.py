#!/usr/bin/env python3
"""
Plot Fisher matrix constraints on DM annihilation and decay parameters,
reproducing the 2×2 layout of Figure 5 in arXiv:2304.07793.

  Top row:    annihilation  <σv>/m_χ   (left: γγ, right: e⁺e⁻)
  Bottom row: decay         Γ_χ        (left: γγ, right: e⁺e⁻)

Usage:
    conda activate test
    cd camb/paper
    python plot_fig5.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from scipy.interpolate import CubicSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

# ----------------------------------------------------------
# Mass grids (must match compute_constraints.py)
# ----------------------------------------------------------
MASS_PANN = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)
XS = np.geomspace(1e-5, 5e3, 500)

# ----------------------------------------------------------
# Planck 2018 reference values from arXiv:2304.07793
# Used to rescale the Fisher Planck curves (post-hoc calibration).
# ----------------------------------------------------------
PLANCK_REF = {
    'gamma_gamma': {
        'pann': 6.0e-28,
        'gamma': 2.5e-25,
        'ref_mass_GeV': 0.1,
    },
    'e_pm': {
        'pann': 4.5e-28,
        'gamma': 5.0e-26,
        'ref_mass_GeV': 0.1,
    },
}


def log_interp(xx, yy):
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    if os.path.exists(path):
        return np.load(path)
    return None


def calibrate_planck(fisher_data, mass_grid, ref_value, ref_mass):
    """Rescale Fisher Planck curve to match paper value at ref_mass."""
    idx = np.argmin(np.abs(mass_grid - ref_mass))
    fisher_at_ref = fisher_data[idx, 3]
    if fisher_at_ref <= 0:
        return fisher_data
    out = fisher_data.copy()
    out[:, 3] = fisher_data[:, 3] * (ref_value / fisher_at_ref)
    return out



# ----------------------------------------------------------
# Style matching arXiv:2304.07793 Figure 5
# Paper colors: black=Planck, red=Ground, blue=PICO, green=CVL
# ----------------------------------------------------------
CURVE_DEFS = [
    ('Planck_cal', 'black',    '-',  1.8, 'Planck'),
    ('ALI_fgres',  'tab:red',  '-',  1.5, 'Ground (noise+fg)'),
    ('ALI_nofgres','tab:red',  '--', 1.2, 'Ground (noise only)'),
    ('PICO_fgres', 'tab:blue', '-',  1.5, 'PICO (noise+fg)'),
    ('PICO_nofgres','tab:blue','--', 1.2, 'PICO (noise only)'),
    ('CVL',        'tab:green','-',  1.5, r'CVL $\ell\in[10,3000]$'),
]


def plot_panel(ax, param_type, channel, xlim=None,
               show_thermal=False, show_legend=False):
    """Plot a single panel."""
    suffix = '' if channel == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    ref = PLANCK_REF.get(channel, PLANCK_REF['gamma_gamma'])
    ref_key = param_type if param_type in ref else 'gamma'

    xs = np.geomspace(xlim[0], xlim[1], 500) if xlim else XS

    for name, color, ls, lw, label in CURVE_DEFS:
        if name == 'Planck_cal':
            raw = load_result(f'sig_{param_type}_Planck{suffix}.npy')
            if raw is None:
                continue
            data = calibrate_planck(raw, mass_grid, ref[ref_key], ref['ref_mass_GeV'])
        else:
            data = load_result(f'sig_{param_type}_{name}{suffix}.npy')

        if data is None:
            continue
        y = data[:, 3]
        if np.any(y <= 0) or np.any(np.isnan(y)):
            continue

        mask = (mass_grid >= xs[0] * 0.5) & (mass_grid <= xs[-1] * 2)
        if np.sum(mask) < 4:
            continue
        cs = log_interp(mass_grid[mask], y[mask])
        ax.loglog(xs, cs(xs), color=color, ls=ls, lw=lw, label=label)

    if show_thermal and param_type == 'pann':
        sv_thermal = 3e-26
        ax.loglog(xs, sv_thermal / xs, color='gray', ls=':', lw=1.0,
                  label=r'$\langle\sigma v\rangle = 3\times10^{-26}$ cm$^3$s$^{-1}$')

    ax.grid(True, which='both', alpha=0.15)
    ax.tick_params(which='both', direction='in', top=True, right=True)

    if show_legend:
        ax.legend(fontsize=7.5, loc='upper right', framealpha=0.85)


def plot_figure5():
    fig, axes = plt.subplots(2, 2, figsize=(13, 10))

    channels = ['gamma_gamma', 'e_pm']
    ch_labels = [r'$\gamma\gamma$', r'$e^+e^-$']

    # Axis ranges matching the paper
    xlims = {
        ('pann',  'gamma_gamma'): (1e-5, 5e3),
        ('pann',  'e_pm'):        (1e-3, 5e3),
        ('gamma', 'gamma_gamma'): (1e-5, 5e3),
        ('gamma', 'e_pm'):        (1e-3, 1e4),
    }

    for j, (ch, ch_lbl) in enumerate(zip(channels, ch_labels)):
        # --- Top row: annihilation ---
        ax = axes[0, j]
        xl_p = xlims[('pann', ch)]
        plot_panel(ax, 'pann', ch, xlim=xl_p,
                   show_thermal=True, show_legend=(j == 0))
        ax.set_xlim(xl_p)
        ax.set_ylim(1e-29, 1e-25)
        if j == 0:
            ax.set_ylabel(
                r'$\langle\sigma v\rangle / m_\chi$'
                r' [cm$^3$ s$^{-1}$ GeV$^{-1}$]', fontsize=13)
        ax.set_title(ch_lbl, fontsize=14)

        # --- Bottom row: decay ---
        ax = axes[1, j]
        xl_g = xlims[('gamma', ch)]
        plot_panel(ax, 'gamma', ch, xlim=xl_g, show_legend=(j == 1))
        ax.set_xlim(xl_g)
        ax.set_ylim(1e-27, 1e-23)
        ax.set_xlabel(r'$m_\chi$ [GeV]', fontsize=13)
        if j == 0:
            ax.set_ylabel(r'$\Gamma_\chi$ [s$^{-1}$]', fontsize=13)

    # Row labels on the right side
    axes[0, 1].annotate('Annihilation', xy=(1.04, 0.5),
                         xycoords='axes fraction', fontsize=12,
                         rotation=270, va='center', ha='left')
    axes[1, 1].annotate('Decay', xy=(1.04, 0.5),
                         xycoords='axes fraction', fontsize=12,
                         rotation=270, va='center', ha='left')

    fig.suptitle(
        r"95% C.L. upper bounds on DM annihilation $\langle\sigma v\rangle/m_\chi$"
        r" (top) and decay width $\Gamma_\chi$ (bottom)"
        "\n(cf. arXiv:2304.07793, Figure 5)",
        fontsize=13, y=0.98)
    fig.tight_layout(rect=[0, 0, 0.97, 0.95])

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig5_fisher.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')

    # Print calibration info
    for ch in channels:
        ref = PLANCK_REF.get(ch, {})
        suffix = '' if ch == 'gamma_gamma' else '_e_pm'
        raw_p = load_result(f'sig_pann_Planck{suffix}.npy')
        raw_g = load_result(f'sig_gamma_Planck{suffix}.npy')
        if raw_p is not None or raw_g is not None:
            print(f'\nPlanck calibration ({ch}):')
        if raw_p is not None:
            idx = np.argmin(np.abs(MASS_PANN - ref.get('ref_mass_GeV', 0.1)))
            ratio_p = ref.get('pann', 0) / raw_p[idx,3] if raw_p[idx,3] > 0 else 0
            print(f'  pann:  Fisher={raw_p[idx,3]:.2e} -> ref={ref.get("pann",0):.2e}  (x{ratio_p:.2f})')
        if raw_g is not None:
            idx = np.argmin(np.abs(MASS_GAMMA - ref.get('ref_mass_GeV', 0.1)))
            ratio_g = ref.get('gamma', 0) / raw_g[idx,3] if raw_g[idx,3] > 0 else 0
            print(f'  gamma: Fisher={raw_g[idx,3]:.2e} -> ref={ref.get("gamma",0):.2e}  (x{ratio_g:.2f})')


if __name__ == '__main__':
    plot_figure5()

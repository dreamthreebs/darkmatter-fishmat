#!/usr/bin/env python3
"""
Plot Fisher matrix constraints on DM annihilation and decay parameters,
reproducing Figure 5 of arXiv:2304.07793 (gamma-gamma channel).

Planck curves use the Fisher shape but are rescaled to match the actual
Planck 2018 results from the paper (post-hoc calibration).

Usage:
    conda activate test
    cd camb/paper
    python plot_fig5.py

Requires: results/ directory populated by compute_constraints.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)


def log_interp(xx, yy):
    """Cubic spline interpolation in log-log space."""
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


# Mass grids (must match compute_constraints.py)
MASS_PANN = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)
XS = np.geomspace(1e-5, 5e3, 500)

# ----------------------------------------------------------
# Planck 2018 reference values from arXiv:2304.07793
# Used to rescale the Fisher Planck curves (post-hoc calibration).
# To switch channels, update these values.
# ----------------------------------------------------------
PLANCK_REF = {
    'gamma_gamma': {
        'pann': 5.6e-28,   # <sv>/m_chi at m~0.1 GeV (paper Sec 4.2)
        'gamma': 3.7e-25,  # Gamma_chi at m~0.1 GeV (~10x PICO, paper conclusion)
        'ref_mass_GeV': 0.1,
    },
    # Placeholder for e+e- channel — fill with paper values when needed
    # 'e_pm': {
    #     'pann': ...,
    #     'gamma': ...,
    #     'ref_mass_GeV': 0.1,
    # },
}


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    if os.path.exists(path):
        return np.load(path)
    return None


def calibrate_planck(fisher_data, mass_grid, ref_value, ref_mass):
    """Rescale Fisher Planck curve so it matches the paper value at ref_mass.
    Preserves the mass-dependent shape from Fisher."""
    idx = np.argmin(np.abs(mass_grid - ref_mass))
    fisher_at_ref = fisher_data[idx, 3]
    if fisher_at_ref <= 0:
        return fisher_data
    scale = ref_value / fisher_at_ref
    out = fisher_data.copy()
    out[:, 3] = fisher_data[:, 3] * scale
    return out


def plot_figure5(channel='gamma_gamma'):
    """Main plotting function — 2 panels (annihilation + decay)."""

    ref = PLANCK_REF[channel]
    fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

    style = {
        'Planck_cal':   {'color': 'black',     'ls': '-',  'lw': 1.8, 'label': 'Planck 2018 (calibrated)'},
        'ALI_fgres':    {'color': 'tab:red',   'ls': '-',  'lw': 1.5, 'label': 'Ground/ALI (noise+fg)'},
        'ALI_nofgres':  {'color': 'tab:red',   'ls': '--', 'lw': 1.2, 'label': 'Ground/ALI (noise only)'},
        'PICO_fgres':   {'color': 'tab:blue',  'ls': '-',  'lw': 1.5, 'label': 'PICO (noise+fg)'},
        'PICO_nofgres': {'color': 'tab:blue',  'ls': '--', 'lw': 1.2, 'label': 'PICO (noise only)'},
        'CVL':          {'color': 'tab:cyan',  'ls': '-',  'lw': 1.5, 'label': r'CVL $\ell\in[10,3000]$'},
    }

    plot_order = ['Planck_cal', 'ALI_fgres', 'ALI_nofgres',
                  'PICO_fgres', 'PICO_nofgres', 'CVL']

    # ----------------------------------------------------------
    # Left panel: annihilation <sigma v>/m_chi
    # ----------------------------------------------------------
    ax = axes[0]
    for name in plot_order:
        if name == 'Planck_cal':
            raw = load_result('sig_pann_Planck.npy')
            if raw is None:
                continue
            data = calibrate_planck(raw, MASS_PANN, ref['pann'], ref['ref_mass_GeV'])
        else:
            data = load_result(f'sig_pann_{name}.npy')
        if data is None:
            print(f'  [skip] sig_pann_{name}.npy not found')
            continue
        y = data[:, 3]
        if np.any(y <= 0) or np.any(np.isnan(y)):
            print(f'  [warn] sig_pann_{name} has invalid values, skipping')
            continue
        cs = log_interp(MASS_PANN, y)
        s = style[name]
        ax.loglog(XS, cs(XS), color=s['color'], ls=s['ls'],
                  lw=s['lw'], label=s['label'])

    sv_thermal = 3e-26
    ax.loglog(XS, sv_thermal / XS, color='gray', ls=':', lw=1.0,
              label=r'$\langle\sigma v\rangle = 3\times10^{-26}$ cm$^3$s$^{-1}$')

    ax.set_xlabel(r'$m_\chi$ [GeV]', fontsize=13)
    ax.set_ylabel(r'$\langle\sigma v\rangle / m_\chi$ [cm$^3$ s$^{-1}$ GeV$^{-1}$]',
                  fontsize=13)
    ax.set_xlim(1e-5, 5e3)
    ax.set_ylim(1e-29, 1e-25)
    ax.set_title(r'Annihilation: $\chi\chi \to \gamma\gamma$', fontsize=14)
    ax.legend(fontsize=8.5, loc='upper right')
    ax.grid(True, which='both', alpha=0.2)

    # ----------------------------------------------------------
    # Right panel: decay Gamma_chi
    # ----------------------------------------------------------
    ax = axes[1]
    for name in plot_order:
        if name == 'Planck_cal':
            raw = load_result('sig_gamma_Planck.npy')
            if raw is None:
                continue
            data = calibrate_planck(raw, MASS_GAMMA, ref['gamma'], ref['ref_mass_GeV'])
        else:
            data = load_result(f'sig_gamma_{name}.npy')
        if data is None:
            print(f'  [skip] sig_gamma_{name}.npy not found')
            continue
        y = data[:, 3]
        if np.any(y <= 0) or np.any(np.isnan(y)):
            print(f'  [warn] sig_gamma_{name} has invalid values, skipping')
            continue
        cs = log_interp(MASS_GAMMA, y)
        s = style[name]
        ax.loglog(XS, cs(XS), color=s['color'], ls=s['ls'],
                  lw=s['lw'], label=s['label'])

    ax.set_xlabel(r'$m_\chi$ [GeV]', fontsize=13)
    ax.set_ylabel(r'$\Gamma_\chi$ [s$^{-1}$]', fontsize=13)
    ax.set_xlim(1e-5, 5e3)
    ax.set_ylim(1e-27, 1e-23)
    ax.set_title(r'Decay: $\chi \to \gamma\gamma$', fontsize=14)
    ax.legend(fontsize=8.5, loc='upper right')
    ax.grid(True, which='both', alpha=0.2)

    # ----------------------------------------------------------
    ch_label = r'$\gamma\gamma$' if channel == 'gamma_gamma' else r'$e^+e^-$'
    fig.suptitle('Fisher Matrix 95% C.L. Upper Bounds on DM Parameters\n'
                 rf'(cf. arXiv:2304.07793, Figure 5, {ch_label} channel)',
                 fontsize=14, y=1.02)
    fig.tight_layout()

    outpath = os.path.join(FIGURE_DIR, 'fig5_fisher.pdf')
    fig.savefig(outpath, dpi=300, bbox_inches='tight')
    print(f'\nFigure saved: {outpath}')

    outpath_png = os.path.join(FIGURE_DIR, 'fig5_fisher.png')
    fig.savefig(outpath_png, dpi=200, bbox_inches='tight')
    print(f'Figure saved: {outpath_png}')

    # Print calibration info
    ref_m = ref['ref_mass_GeV']
    raw_pann = load_result('sig_pann_Planck.npy')
    raw_gamma = load_result('sig_gamma_Planck.npy')
    if raw_pann is not None:
        idx = np.argmin(np.abs(MASS_PANN - ref_m))
        print(f'\nPlanck calibration (m_ref={ref_m} GeV):')
        print(f'  pann:  Fisher={raw_pann[idx,3]:.2e} -> paper={ref["pann"]:.2e}  (x{ref["pann"]/raw_pann[idx,3]:.2f})')
    if raw_gamma is not None:
        idx = np.argmin(np.abs(MASS_GAMMA - ref_m))
        print(f'  gamma: Fisher={raw_gamma[idx,3]:.2e} -> paper={ref["gamma"]:.2e}  (x{ref["gamma"]/raw_gamma[idx,3]:.2f})')


if __name__ == '__main__':
    plot_figure5()

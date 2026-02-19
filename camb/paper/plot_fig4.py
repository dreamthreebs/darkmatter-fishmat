#!/usr/bin/env python3
"""
Plot Figure 4 of arXiv:2304.07793:
  TT vs EE vs TT+TE+EE constraints for ground observation (ALI), γγ channel.

  Left:  annihilation  σv/m_χ
  Right: decay         Γ_χ

Usage:
    cd camb/paper
    python plot_fig4.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from scipy.interpolate import CubicSpline, UnivariateSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

MASS_PANN = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)


def log_interp(xx, yy, smooth=0):
    if smooth > 0:
        sp = UnivariateSpline(np.log10(xx), np.log10(yy), s=smooth)
        return lambda zz: np.power(10, sp(np.log10(zz)))
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    if os.path.exists(path):
        return np.load(path)
    return None


CURVE_DEFS = [
    # (col_index, color, ls, lw, label)  — col: 0=TT, 1=EE, 3=TT+TE+EE
    (0, 'black',    '-',  1.8, 'TT'),
    (1, 'tab:red',  '-',  1.8, 'EE'),
    (3, 'tab:blue', '-',  1.8, 'TT+TE+EE'),
]


def plot_figure4():
    fig, axes = plt.subplots(1, 2, figsize=(13, 5))

    panels = [
        ('pann',  MASS_PANN,  'sig_pann_ALI_fgres.npy',
         r'$\langle\sigma v\rangle / m_\chi$'
         r' [cm$^3$ s$^{-1}$ GeV$^{-1}$]',
         (1e-29, 1e-24)),
        ('gamma', MASS_GAMMA, 'sig_gamma_ALI_fgres.npy',
         r'$\Gamma_\chi$ [s$^{-1}$]',
         (1e-27, 1e-22)),
    ]

    for ax, (ptype, mass_grid, fname, ylabel, ylim) in zip(axes, panels):
        data = load_result(fname)
        if data is None:
            ax.text(0.5, 0.5, f'{fname} not found', transform=ax.transAxes,
                    ha='center')
            continue

        xs = np.geomspace(mass_grid[0], mass_grid[-1], 500)

        for col, color, ls, lw, label in CURVE_DEFS:
            y = data[:, col]
            if np.any(y <= 0) or np.any(np.isnan(y)):
                continue
            sm = 0.15
            cs = log_interp(mass_grid, y, smooth=sm)
            ax.loglog(xs, cs(xs), color=color, ls=ls, lw=lw, label=label)

        ax.set_xlim(mass_grid[0], mass_grid[-1])
        ax.set_ylim(ylim)
        ax.set_xlabel(r'$m_\chi$ [GeV]', fontsize=13)
        ax.set_ylabel(ylabel, fontsize=13)
        ax.grid(True, which='both', alpha=0.15)
        ax.tick_params(which='both', direction='in', top=True, right=True)
        ax.legend(fontsize=11, loc='upper right', framealpha=0.85)

    axes[0].set_title('Annihilation', fontsize=14)
    axes[1].set_title('Decay', fontsize=14)

    fig.suptitle(
        r"95% C.L. upper bounds — Ground observation, $\gamma\gamma$ channel"
        "\n(NILC V1, cf. arXiv:2304.07793 Figure 4)",
        fontsize=13, y=1.02)
    fig.tight_layout()

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig4_TT_EE_combined.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure4()

#!/usr/bin/env python3
"""
Figure 5: Fisher matrix sensitivity on DM annihilation/decay parameters
(γγ channel only, NO Planck line).

  Top row:    Ground  (left: annihilation, right: decay)
  Bottom row: PICO    (left: annihilation, right: decay)

Each panel shows noise+fg, noise-only, and CVL curves.

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
from scipy.interpolate import CubicSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

MASS_PANN = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)

CHANNEL = 'gamma_gamma'

EXPERIMENTS = [
    {
        'label': 'Ground',
        'curves': [
            ('ALI_fgres',   'tab:red',  '-',  1.5, 'noise + foreground'),
            ('ALI_nofgres', 'tab:red',  '--', 1.2, 'noise only'),
            ('CVL',         'black',    '-',  1.5, r'CVL $\ell\in[10,3000]$'),
        ],
    },
    {
        'label': 'PICO',
        'curves': [
            ('PICO_fgres',   'tab:blue', '-',  1.5, 'noise + foreground'),
            ('PICO_nofgres', 'tab:blue', '--', 1.2, 'noise only'),
            ('CVL',          'black',    '-',  1.5, r'CVL $\ell\in[10,3000]$'),
        ],
    },
]


def log_interp(xx, yy):
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    if os.path.exists(path):
        return np.load(path)
    return None


def plot_panel(ax, param_type, curves, xlim=None, show_legend=False):
    suffix = '' if CHANNEL == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500) if xlim else np.geomspace(1e-5, 5e3, 500)

    for name, color, ls, lw, label in curves:
        data = load_result(f'sig_{param_type}_{name}{suffix}.npy')
        if data is None:
            continue
        y = data[:, 3]
        if np.any(y <= 0) or np.any(np.isnan(y)):
            continue
        mask = (mass_grid >= xs[0] * 0.5) & (mass_grid <= xs[-1] * 2)
        if np.sum(mask) < 4:
            continue
        interp = log_interp(mass_grid[mask], y[mask])
        ax.loglog(xs, interp(xs), color=color, ls=ls, lw=lw, label=label)

    import matplotlib.ticker as mticker
    ax.xaxis.set_major_locator(mticker.LogLocator(base=10, numticks=20))
    ax.xaxis.set_major_formatter(mticker.LogFormatterSciNotation())
    ax.xaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    ax.xaxis.set_minor_formatter(mticker.NullFormatter())
    ax.yaxis.set_minor_locator(mticker.LogLocator(base=10, subs=np.arange(2, 10), numticks=100))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    ax.grid(True, which='major', alpha=0.2)
    ax.grid(True, which='minor', alpha=0.08)
    ax.tick_params(which='major', direction='in', top=True, right=True)
    ax.tick_params(which='minor', direction='in', top=True, right=True, length=3)
    if show_legend:
        ax.legend(fontsize=7.5, loc='upper right', framealpha=0)


def plot_figure5():
    fig, axes = plt.subplots(2, 2, figsize=(9, 7))

    col_titles = [r'$\chi\chi\to\gamma\gamma$', r'$\chi\to\gamma\gamma$']
    param_types = ['pann', 'gamma']
    xlims = {'pann': (1e-5, 5e3), 'gamma': (1e-5, 5e3)}
    ylims = {'pann': (5e-29, 5e-27), 'gamma': (3e-27, 2e-24)}
    ylabels = {
        'pann':  r'$\langle\sigma v\rangle / m_\chi$'
                 r' [cm$^3$ s$^{-1}$ GeV$^{-1}$]',
        'gamma': r'$\Gamma_\chi$ [s$^{-1}$]',
    }

    for row, exp in enumerate(EXPERIMENTS):
        for col, pt in enumerate(param_types):
            ax = axes[row, col]
            xl = xlims[pt]
            plot_panel(ax, pt, exp['curves'], xlim=xl,
                       show_legend=(col == 0))
            ax.set_xlim(xl)
            ax.set_ylim(ylims[pt])
            if col == 0:
                ax.set_ylabel(ylabels[pt], fontsize=12)
            if row == 1:
                ax.set_xlabel(r'$m_\chi$ [GeV]', fontsize=13)
            if row == 0:
                ax.set_title(col_titles[col], fontsize=14)

    for row, exp in enumerate(EXPERIMENTS):
        axes[row, 1].annotate(exp['label'], xy=(1.04, 0.5),
                              xycoords='axes fraction', fontsize=12,
                              rotation=270, va='center', ha='left')

    fig.tight_layout(rect=[0, 0, 0.97, 1])

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig5_fisher.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure5()

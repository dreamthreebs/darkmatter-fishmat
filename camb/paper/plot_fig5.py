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
from plot_style import (
    FONT_SIZE,
    LEGEND_SIZE,
    LINE_WIDTH,
    LINE_WIDTH_THIN,
    TITLE_SIZE,
    apply_legend_style,
    apply_result_axis_style,
    configure_matplotlib,
)

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
            ('ALI_fgres',   'tab:red',  '-',  LINE_WIDTH, 'noise + foreground'),
            ('ALI_nofgres', 'tab:red',  '--', LINE_WIDTH_THIN, 'noise only'),
            ('CVL',         'black',    '-',  LINE_WIDTH, r'CVL $\ell\in[10,3000]$'),
        ],
    },
    {
        'label': 'PICO',
        'curves': [
            ('PICO_fgres',   'tab:blue', '-',  LINE_WIDTH, 'noise + foreground'),
            ('PICO_nofgres', 'tab:blue', '--', LINE_WIDTH_THIN, 'noise only'),
            ('CVL',          'black',    '-',  LINE_WIDTH, r'CVL $\ell\in[10,3000]$'),
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


def auto_log_ylim(y_arrays, pad_fraction=0.25, min_span_decades=1.0):
    if not y_arrays:
        return None
    y = np.concatenate(y_arrays)
    y = y[np.isfinite(y) & (y > 0)]
    if y.size == 0:
        return None

    log_y = np.log10(y)
    lo = np.min(log_y)
    hi = np.max(log_y)
    center = 0.5 * (lo + hi)
    span = max(hi - lo, min_span_decades)
    span *= 1.0 + 2.0 * pad_fraction
    return 10.0 ** (center - 0.5 * span), 10.0 ** (center + 0.5 * span)


def plot_panel(ax, param_type, curves, xlim=None, show_legend=False):
    suffix = '' if CHANNEL == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500) if xlim else np.geomspace(1e-5, 5e3, 500)
    plotted_y = []
    if xlim:
        ax.set_xlim(xlim)

    for name, color, ls, lw, label in curves:
        if name == 'CVL':
            continue
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
        yy = interp(xs)
        ax.plot(xs, yy, color=color, ls=ls, lw=lw, label=label)
        plotted_y.append(yy)

    apply_result_axis_style(ax)
    if show_legend:
        apply_legend_style(ax, fontsize=LEGEND_SIZE, loc='lower right')
    return auto_log_ylim(plotted_y)


def plot_figure5():
    configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(24, 24))

    col_titles = [r'$\chi\chi \rightarrow \gamma\gamma$', r'$\chi \rightarrow \gamma\gamma$']
    param_types = ['pann', 'gamma']
    xlims = {'pann': (1e-5, 5e3), 'gamma': (1e-5, 5e3)}
    ylabels = {
        'pann':  r'$\langle\sigma v\rangle/m_\chi\quad[cm^3s^{-1}GeV^{-1}]$',
        'gamma': r'$\Gamma_\chi[s^{-1}]$',
    }

    for row, exp in enumerate(EXPERIMENTS):
        for col, pt in enumerate(param_types):
            ax = axes[row, col]
            xl = xlims[pt]
            ylim = plot_panel(ax, pt, exp['curves'], xlim=xl,
                              show_legend=(col == 1))
            ax.set_xlim(xl)
            if ylim is not None:
                ax.set_ylim(ylim)
            ax.set_ylabel(ylabels[pt], fontsize=FONT_SIZE)
            ax.set_xlabel(r'$m_\chi[GeV]$', fontsize=FONT_SIZE)
            if row == 0:
                ax.set_title(col_titles[col], fontsize=TITLE_SIZE, pad=10)

    for row, exp in enumerate(EXPERIMENTS):
        axes[row, 1].annotate(exp['label'], xy=(1.04, 0.5),
                              xycoords='axes fraction', fontsize=FONT_SIZE,
                              rotation=270, va='center', ha='left')

    fig.tight_layout(rect=[0, 0, 0.97, 0.96])

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig5_fisher.{fmt}')
        fig.savefig(path, dpi=dpi)
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure5()

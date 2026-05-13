#!/usr/bin/env python3
"""
Figure 6 (v2 with ratio panels).

Main panels (top):    sigma upper limits from Planck (calibrated), Ground (ALI),
                      PICO, CVL.
Ratio panels (below): improvement factor   sigma_Planck / sigma_{experiment} ,
                      i.e. "how many times tighter the bound becomes when going
                      from Planck to Ground/PICO/CVL".  Shows experiment hierarchy
                      and the remaining headroom up to CVL.

2x2 layout (each cell is main + ratio):
  rows:    pann (annihilation, <sigma v>/m_chi)   |   gamma (decay, Gamma_chi)
  cols:    gamma-gamma                            |   e+ e-
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import matplotlib.ticker as mticker
from scipy.interpolate import CubicSpline

from plot_style_v2 import (
    FONT_SIZE,
    LEGEND_SIZE,
    LINE_WIDTH,
    TITLE_SIZE,
    C_PLANCK,
    C_GROUND,
    C_PICO,
    C_CVL,
    apply_axis_style_v2,
    apply_ratio_axis_style_v2,
    apply_joined_pair_x_style,
    apply_legend_style_v2,
    configure_matplotlib_v2,
)

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures_v2')
os.makedirs(FIGURE_DIR, exist_ok=True)

MASS_PANN = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)

PLANCK_REF = {
    'gamma_gamma': {'pann': 6.0e-28, 'gamma': 2.5e-25, 'ref_mass_GeV': 0.1},
    'e_pm':        {'pann': 4.5e-28, 'gamma': 5.0e-26, 'ref_mass_GeV': 0.1},
}


def log_interp(xx, yy):
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    return np.load(path) if os.path.exists(path) else None


def calibrate_planck(fisher_data, mass_grid, ref_value, ref_mass):
    idx = np.argmin(np.abs(mass_grid - ref_mass))
    fisher_at_ref = fisher_data[idx, 3]
    if fisher_at_ref <= 0:
        return fisher_data
    out = fisher_data.copy()
    out[:, 3] = fisher_data[:, 3] * (ref_value / fisher_at_ref)
    return out


CURVE_DEFS = [
    ('Planck_cal', C_PLANCK, '-', LINE_WIDTH, 'Planck 2018'),
    ('ALI_fgres',  C_GROUND, '-', LINE_WIDTH, 'Ground (ALI)'),
    ('PICO_fgres', C_PICO,   '-', LINE_WIDTH, 'PICO'),
    ('CVL',        C_CVL,    '-', LINE_WIDTH, r'CVL $\ell\in[10,3000]$'),
]


def get_curve_y(name, param_type, channel, mass_grid):
    suffix = '' if channel == 'gamma_gamma' else '_e_pm'
    ref = PLANCK_REF.get(channel, PLANCK_REF['gamma_gamma'])
    ref_key = param_type if param_type in ref else 'gamma'
    if name == 'Planck_cal':
        raw = load_result(f'sig_{param_type}_Planck{suffix}.npy')
        if raw is None:
            return None
        return calibrate_planck(raw, mass_grid, ref[ref_key], ref['ref_mass_GeV'])[:, 3]
    data = load_result(f'sig_{param_type}_{name}{suffix}.npy')
    if data is None:
        return None
    return data[:, 3]


def plot_main(ax, param_type, channel, xlim, show_legend=False):
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500)

    for name, color, ls, lw, label in CURVE_DEFS:
        y = get_curve_y(name, param_type, channel, mass_grid)
        if y is None:
            continue
        if np.any(y <= 0) or np.any(np.isnan(y)):
            continue
        mask = (mass_grid >= xs[0] * 0.5) & (mass_grid <= xs[-1] * 2)
        if np.sum(mask) < 4:
            continue
        interp = log_interp(mass_grid[mask], y[mask])
        ax.plot(xs, interp(xs), color=color, ls=ls, lw=lw, label=label,
                solid_capstyle='round', dash_capstyle='round')

    apply_axis_style_v2(ax)
    if show_legend:
        apply_legend_style_v2(ax, fontsize=LEGEND_SIZE, loc='lower right')


def plot_ratio(ax, param_type, channel, xlim, show_legend=False):
    """Each non-CVL curve divided by CVL: how many times *weaker* than the
    CVL bound it is.  Larger ratio = larger gap from the ideal limit.
    """
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500)

    y_cvl = get_curve_y('CVL', param_type, channel, mass_grid)
    if y_cvl is None:
        return

    for name, color, _ls, _lw, label in CURVE_DEFS:
        if name == 'CVL':
            continue
        y = get_curve_y(name, param_type, channel, mass_grid)
        if y is None:
            continue
        good = (y_cvl > 0) & (y > 0)
        if np.sum(good) < 4:
            continue
        ratio = y[good] / y_cvl[good]
        cs = log_interp(mass_grid[good], ratio)
        ax.plot(xs, cs(xs), color=color, lw=LINE_WIDTH, label=label,
                solid_capstyle='round')

    ax.set_xlim(xlim)
    if param_type == 'pann':
        ax.set_ylim(0.7, 10)
        yticks = [1, 2, 5, 10]
    else:
        ax.set_ylim(0.7, 30)
        yticks = [1, 2, 5, 10, 30]
    apply_ratio_axis_style_v2(
        ax,
        ylabel='Ratio to CVL',
        y_is_log=True,
        yticks=yticks,
        hline=1.0,
    )
    ax.yaxis.set_major_formatter(mticker.FuncFormatter(lambda v, _: f'{v:g}'))
    ax.yaxis.set_minor_formatter(mticker.NullFormatter())
    if show_legend:
        apply_legend_style_v2(ax, fontsize=LEGEND_SIZE - 6,
                              loc='upper right', ncol=1)


def plot_figure6_v2():
    configure_matplotlib_v2()

    fig = plt.figure(figsize=(22, 18), constrained_layout=False)
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        hspace=0.22,
        left=0.075, right=0.97, top=0.95, bottom=0.06,
    )
    inner_blocks = [
        gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=outer[i],
            height_ratios=[3.0, 1.0], width_ratios=[1, 1],
            hspace=0.0, wspace=0.22,
        )
        for i in range(2)
    ]

    channels = ['gamma_gamma', 'e_pm']
    ch_titles = {
        ('pann', 'gamma_gamma'): r'$\chi\chi \to \gamma\gamma$',
        ('pann', 'e_pm'):        r'$\chi\chi \to e^- e^+$',
        ('gamma', 'gamma_gamma'): r'$\chi \to \gamma\gamma$',
        ('gamma', 'e_pm'):        r'$\chi \to e^- e^+$',
    }
    xlims = {
        ('pann',  'gamma_gamma'): (1e-5, 5e3),
        ('pann',  'e_pm'):        (1e-3, 5e3),
        ('gamma', 'gamma_gamma'): (1e-5, 5e3),
        ('gamma', 'e_pm'):        (1e-3, 1e4),
    }
    main_ylims = {'pann': (1e-29, 1e-25), 'gamma': (1e-27, 1e-23)}
    main_ylabels = {
        'pann':  r'$\langle\sigma v\rangle / m_\chi\ \ [\mathrm{cm^{3}\,s^{-1}\,GeV^{-1}}]$',
        'gamma': r'$\Gamma_\chi\ \ [\mathrm{s^{-1}}]$',
    }

    main_axes = [[None, None], [None, None]]

    for row_block, pt in enumerate(['pann', 'gamma']):
        inner = inner_blocks[row_block]
        for col, ch in enumerate(channels):
            ax_main = fig.add_subplot(inner[0, col])
            ax_ratio = fig.add_subplot(inner[1, col], sharex=ax_main)

            xl = xlims[(pt, ch)]
            plot_main(ax_main, pt, ch, xlim=xl,
                      show_legend=(row_block == 0 and col == 1))
            ax_main.set_xlim(xl)
            ax_main.set_ylim(main_ylims[pt])
            ax_main.set_ylabel(main_ylabels[pt], fontsize=FONT_SIZE)
            ax_main.set_title(ch_titles[(pt, ch)], fontsize=TITLE_SIZE, pad=12)

            plot_ratio(ax_ratio, pt, ch, xlim=xl, show_legend=False)
            ax_ratio.set_xlabel(r'$m_\chi\ \ [\mathrm{GeV}]$', fontsize=FONT_SIZE)

            apply_joined_pair_x_style(ax_main, ax_ratio)

            main_axes[row_block][col] = (ax_main, ax_ratio)

    panel_tags = [['(a)', '(b)'], ['(c)', '(d)']]
    for r in range(2):
        for c in range(2):
            ax = main_axes[r][c][0]
            ax.text(0.025, 0.97, panel_tags[r][c],
                    transform=ax.transAxes, ha='left', va='top',
                    fontsize=TITLE_SIZE, fontweight='bold', color='#222222')

    for fmt, dpi in [('pdf', 300), ('png', 220)]:
        path = os.path.join(FIGURE_DIR, f'fig6_fisher_v2.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure6_v2()

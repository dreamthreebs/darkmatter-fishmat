#!/usr/bin/env python3
"""
Figure 5 (v2 with degradation panels).

Top of each column: noise+foreground (solid) vs noise only (dashed) vs CVL.
Degradation panel:  (sigma_with_fg / sigma_noise_only - 1) * 100 %,
                    i.e. how much the foreground residual loosens the bound.

Rows: Ground (ALI) above, PICO below.
Columns: annihilation <sigma v>/m_chi  &  decay  Gamma_chi  (gamma-gamma).
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
    LINE_WIDTH_THIN,
    TITLE_SIZE,
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

CHANNEL = 'gamma_gamma'

EXPERIMENTS = [
    {
        'label': 'Ground (ALI)',
        'tag_color': C_GROUND,
        'curves': [
            ('ALI_fgres',   C_GROUND, '-',  LINE_WIDTH,      'noise + foreground'),
            ('ALI_nofgres', C_GROUND, '--', LINE_WIDTH_THIN, 'noise only'),
            ('CVL',         C_CVL,    '-',  LINE_WIDTH,      r'CVL $\ell\in[10,3000]$'),
        ],
    },
    {
        'label': 'PICO',
        'tag_color': C_PICO,
        'curves': [
            ('PICO_fgres',   C_PICO,  '-',  LINE_WIDTH,      'noise + foreground'),
            ('PICO_nofgres', C_PICO,  '--', LINE_WIDTH_THIN, 'noise only'),
            ('CVL',          C_CVL,   '-',  LINE_WIDTH,      r'CVL $\ell\in[10,3000]$'),
        ],
    },
]


def log_interp(xx, yy):
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    return np.load(path) if os.path.exists(path) else None


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


def plot_main(ax, param_type, curves, xlim, show_legend=False):
    suffix = '' if CHANNEL == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500)
    plotted_y = []
    ax.set_xlim(xlim)

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
        yy = interp(xs)
        ax.plot(xs, yy, color=color, ls=ls, lw=lw, label=label,
                solid_capstyle='round', dash_capstyle='round')
        plotted_y.append(yy)

    apply_axis_style_v2(ax)
    if show_legend:
        apply_legend_style_v2(ax, fontsize=LEGEND_SIZE, loc='lower right')
    return auto_log_ylim(plotted_y)


def plot_ratio(ax, param_type, curves, xlim, show_legend=False):
    """Foreground-induced constraint degradation in percent:

        delta_sigma = (sigma_with_fg / sigma_noise_only - 1) * 100 %.

    One curve per experiment (color inherited from the main panel).
    """
    suffix = '' if CHANNEL == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    xs = np.geomspace(xlim[0], xlim[1], 500)

    fg_name = nofg_name = None
    exp_color = None
    for name, color, _ls, _lw, _label in curves:
        if name.endswith('_fgres'):
            fg_name = name
            exp_color = color
        elif name.endswith('_nofgres'):
            nofg_name = name
    if fg_name is None or nofg_name is None or exp_color is None:
        return

    d_fg = load_result(f'sig_{param_type}_{fg_name}{suffix}.npy')
    d_no = load_result(f'sig_{param_type}_{nofg_name}{suffix}.npy')
    if d_fg is None or d_no is None:
        return

    good = (d_fg[:, 3] > 0) & (d_no[:, 3] > 0)
    if np.sum(good) < 4:
        return
    deg_pct = (d_fg[good, 3] / d_no[good, 3] - 1.0) * 100.0
    cs = CubicSpline(np.log10(mass_grid[good]), deg_pct)
    ax.plot(xs, cs(np.log10(xs)), color=exp_color, ls='-', lw=LINE_WIDTH,
            solid_capstyle='round')

    ax.set_xlim(xlim)
    ax.set_ylim(0, 35)
    apply_ratio_axis_style_v2(
        ax,
        ylabel=r'$\dfrac{\mathrm{with\ FG}}{\mathrm{no\ FG}}-1\ \ [\%]$',
        y_is_log=False,
        hline=0.0,
    )
    ax.yaxis.set_major_locator(mticker.MultipleLocator(10))
    ax.yaxis.set_minor_locator(mticker.MultipleLocator(5))
    ax.yaxis.set_major_formatter(
        mticker.FuncFormatter(lambda v, _: f'{int(round(v))}')
    )
    ax.yaxis.label.set_fontsize(FONT_SIZE - 12)


def plot_figure5_v2():
    configure_matplotlib_v2()

    fig = plt.figure(figsize=(22, 18), constrained_layout=False)
    outer = gridspec.GridSpec(
        2, 1, figure=fig,
        hspace=0.22,
        left=0.075, right=0.96, top=0.95, bottom=0.06,
    )
    inner_blocks = [
        gridspec.GridSpecFromSubplotSpec(
            2, 2, subplot_spec=outer[i],
            height_ratios=[3.0, 1.35], width_ratios=[1, 1],
            hspace=0.0, wspace=0.20,
        )
        for i in range(2)
    ]

    col_titles = [r'$\chi\chi \to \gamma\gamma$', r'$\chi \to \gamma\gamma$']
    param_types = ['pann', 'gamma']
    xlims = {'pann': (1e-5, 5e3), 'gamma': (1e-5, 5e3)}
    ylabels = {
        'pann':  r'$\langle\sigma v\rangle / m_\chi\ \ [\mathrm{cm^{3}\,s^{-1}\,GeV^{-1}}]$',
        'gamma': r'$\Gamma_\chi\ \ [\mathrm{s^{-1}}]$',
    }

    main_axes_all = []

    for row_block, exp in enumerate(EXPERIMENTS):
        inner = inner_blocks[row_block]
        main_axes_row = []
        for col, pt in enumerate(param_types):
            ax_main = fig.add_subplot(inner[0, col])
            ax_ratio = fig.add_subplot(inner[1, col], sharex=ax_main)

            xl = xlims[pt]
            ylim = plot_main(
                ax_main, pt, exp['curves'], xlim=xl,
                show_legend=(col == 1 and row_block == 1),
            )
            ax_main.set_xlim(xl)
            if ylim is not None:
                ax_main.set_ylim(ylim)
            ax_main.set_ylabel(ylabels[pt], fontsize=FONT_SIZE)
            ax_main.set_title(col_titles[col], fontsize=TITLE_SIZE, pad=12)

            plot_ratio(ax_ratio, pt, exp['curves'], xlim=xl)
            ax_ratio.set_xlabel(r'$m_\chi\ \ [\mathrm{GeV}]$', fontsize=FONT_SIZE)

            apply_joined_pair_x_style(ax_main, ax_ratio)

            main_axes_row.append((ax_main, ax_ratio))

        # Per-row experiment label on right
        ax_right = main_axes_row[1][0]
        ax_right.annotate(
            exp['label'], xy=(1.04, 0.5), xycoords='axes fraction',
            fontsize=FONT_SIZE, rotation=270, va='center', ha='left',
            fontweight='bold', color=exp['tag_color'],
        )

        main_axes_all.append(main_axes_row)

    panel_tags = [['(a)', '(b)'], ['(c)', '(d)']]
    for row_block in range(2):
        for col in range(2):
            ax = main_axes_all[row_block][col][0]
            ax.text(0.025, 0.97, panel_tags[row_block][col],
                    transform=ax.transAxes, ha='left', va='top',
                    fontsize=TITLE_SIZE, fontweight='bold', color='#222222')

    for fmt, dpi in [('pdf', 300), ('png', 220)]:
        path = os.path.join(FIGURE_DIR, f'fig5_fisher_v2.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure5_v2()

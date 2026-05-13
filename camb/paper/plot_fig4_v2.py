#!/usr/bin/env python3
"""
Figure 4 (v2 with ratio panels).

Main panels  (top):   TT / EE / TT+TE+EE  2sigma upper limits for ALI ground
                      observation, gamma-gamma channel.
Ratio panels (below): TT / EE   and   (TT+TE+EE) / EE   --- shows that EE
                      tightens the bound by ~ an order of magnitude over TT,
                      and that adding TT/TE to EE only gives a ~ factor-2
                      improvement, i.e. the constraint is polarization-driven.

Left column:  annihilation  <sigma v>/m_chi
Right column: decay          Gamma_chi
"""

import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
from scipy.interpolate import CubicSpline, UnivariateSpline

from plot_style_v2 import (
    FONT_SIZE,
    LEGEND_SIZE,
    LINE_WIDTH,
    TITLE_SIZE,
    C_TT,
    C_EE,
    C_COMBINED,
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


def log_interp(xx, yy, smooth=0):
    if smooth > 0:
        sp = UnivariateSpline(np.log10(xx), np.log10(yy), s=smooth)
        return lambda zz: np.power(10, sp(np.log10(zz)))
    cs = CubicSpline(np.log10(xx), np.log10(yy))
    return lambda zz: np.power(10, cs(np.log10(zz)))


def load_result(name):
    path = os.path.join(RESULTS_DIR, name)
    return np.load(path) if os.path.exists(path) else None


CURVE_DEFS = [
    (0, C_TT,       '-', LINE_WIDTH, 'TT'),
    (1, C_EE,       '-', LINE_WIDTH, 'EE'),
    (3, C_COMBINED, '-', LINE_WIDTH, 'TT + TE + EE'),
]


def plot_figure4_v2():
    configure_matplotlib_v2()

    fig = plt.figure(figsize=(22, 11), constrained_layout=False)
    gs = gridspec.GridSpec(
        2, 2, figure=fig,
        height_ratios=[3.0, 1.2],
        width_ratios=[1, 1],
        hspace=0.0, wspace=0.20,
        left=0.075, right=0.97, top=0.93, bottom=0.09,
    )

    panels = [
        ('pann',  MASS_PANN,  'sig_pann_ALI_fgres.npy',
         r'$\langle\sigma v\rangle / m_\chi\ \ [\mathrm{cm^{3}\,s^{-1}\,GeV^{-1}}]$',
         (1e-28, 1e-24),
         r'$\chi\chi \to \gamma\gamma$  (Ground / ALI)'),
        ('gamma', MASS_GAMMA, 'sig_gamma_ALI_fgres.npy',
         r'$\Gamma_\chi\ \ [\mathrm{s^{-1}}]$',
         (1e-26, 1e-22),
         r'$\chi \to \gamma\gamma$  (Ground / ALI)'),
    ]

    main_axes = []
    ratio_axes = []

    for i, (ptype, mass_grid, fname, ylabel, ylim, title) in enumerate(panels):
        ax_main = fig.add_subplot(gs[0, i])
        ax_ratio = fig.add_subplot(gs[1, i], sharex=ax_main)

        data = load_result(fname)
        if data is None:
            ax_main.text(0.5, 0.5, f'{fname} not found',
                         transform=ax_main.transAxes, ha='center')
            main_axes.append(ax_main)
            ratio_axes.append(ax_ratio)
            continue

        xs = np.geomspace(mass_grid[0], mass_grid[-1], 500)

        # Main panel: absolute upper limits
        for col, color, ls, lw, label in CURVE_DEFS:
            y = data[:, col]
            if np.any(y <= 0) or np.any(np.isnan(y)):
                continue
            cs = log_interp(mass_grid, y, smooth=0.15)
            ax_main.plot(xs, cs(xs), color=color, ls=ls, lw=lw, label=label,
                         solid_capstyle='round', dash_capstyle='round')

        ax_main.set_xlim(1e-5, 5e3)
        ax_main.set_ylim(ylim)
        ax_main.set_ylabel(ylabel, fontsize=FONT_SIZE)
        ax_main.set_title(title, fontsize=TITLE_SIZE, pad=12)
        apply_axis_style_v2(ax_main)
        if i == 1:
            apply_legend_style_v2(ax_main, fontsize=LEGEND_SIZE, loc='lower right')

        # Ratio panel: each main curve divided by the strongest curve
        # (TT+TE+EE).  Reads as "factor weaker than the combined bound".
        tt = data[:, 0]; ee = data[:, 1]; comb = data[:, 3]
        good = (tt > 0) & (ee > 0) & (comb > 0)
        m_good = mass_grid[good]

        r_tt = log_interp(m_good, tt[good] / comb[good], smooth=0.05)
        r_ee = log_interp(m_good, ee[good] / comb[good], smooth=0.05)
        ax_ratio.plot(xs, r_tt(xs), color=C_TT, ls='-', lw=LINE_WIDTH,
                      label='TT', solid_capstyle='round')
        ax_ratio.plot(xs, r_ee(xs), color=C_EE, ls='-', lw=LINE_WIDTH,
                      label='EE', solid_capstyle='round')

        ax_ratio.set_xlim(1e-5, 5e3)
        ax_ratio.set_ylim(0.5, 100)
        ax_ratio.set_xlabel(r'$m_\chi\ \ [\mathrm{GeV}]$', fontsize=FONT_SIZE)
        apply_ratio_axis_style_v2(
            ax_ratio,
            ylabel='Ratio to TT+TE+EE',
            y_is_log=True,
            yticks=[1, 3, 10, 30, 100],
            hline=1.0,
        )
        ax_ratio.yaxis.set_major_formatter(
            matplotlib.ticker.FuncFormatter(lambda v, _: f'{v:g}')
        )
        ax_ratio.yaxis.set_minor_formatter(matplotlib.ticker.NullFormatter())

        apply_joined_pair_x_style(ax_main, ax_ratio)

        main_axes.append(ax_main)
        ratio_axes.append(ax_ratio)

    for ax, tag in zip(main_axes, ['(a)', '(b)']):
        ax.text(0.02, 0.97, tag, transform=ax.transAxes,
                ha='left', va='top', fontsize=TITLE_SIZE,
                fontweight='bold', color='#222222')

    for fmt, dpi in [('pdf', 300), ('png', 220)]:
        path = os.path.join(FIGURE_DIR, f'fig4_TT_EE_combined_v2.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure4_v2()

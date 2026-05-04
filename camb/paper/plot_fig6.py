#!/usr/bin/env python3
"""
Figure 6: Fisher matrix constraints compared with Planck 2018.

Same layout as Figure 5 (2×2: annihilation/decay × γγ/e⁺e⁻),
but adds the calibrated Planck line for comparison.

Usage:
    conda activate test
    cd camb/paper
    python plot_fig6.py
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
    idx = np.argmin(np.abs(mass_grid - ref_mass))
    fisher_at_ref = fisher_data[idx, 3]
    if fisher_at_ref <= 0:
        return fisher_data
    out = fisher_data.copy()
    out[:, 3] = fisher_data[:, 3] * (ref_value / fisher_at_ref)
    return out


CURVE_DEFS = [
    # ('Planck_cal', 'black',     '-', LINE_WIDTH, 'Planck'),
    ('ALI_fgres',  'tab:red',   '-', LINE_WIDTH, 'Ground observation'),
    ('PICO_fgres', 'tab:green', '-', LINE_WIDTH, 'PICO'),
    ('CVL',        'tab:blue',  '-', LINE_WIDTH, r'CVL $\ell\in[10,3000]$'),
]


def plot_panel(ax, param_type, channel, xlim=None, show_legend=False):
    suffix = '' if channel == 'gamma_gamma' else '_e_pm'
    mass_grid = MASS_PANN if param_type == 'pann' else MASS_GAMMA
    ref = PLANCK_REF.get(channel, PLANCK_REF['gamma_gamma'])
    ref_key = param_type if param_type in ref else 'gamma'
    xs = np.geomspace(xlim[0], xlim[1], 500) if xlim else np.geomspace(1e-5, 5e3, 500)

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
        interp = log_interp(mass_grid[mask], y[mask])
        ax.plot(xs, interp(xs), color=color, ls=ls, lw=lw, label=label)

    apply_result_axis_style(ax)
    if show_legend:
        apply_legend_style(ax, fontsize=LEGEND_SIZE, loc='lower right')


def plot_figure6():
    configure_matplotlib()
    fig, axes = plt.subplots(2, 2, figsize=(24, 24))

    channels = ['gamma_gamma', 'e_pm']
    ch_titles = {
        ('pann', 'gamma_gamma'): r'$\chi\chi \rightarrow \gamma\gamma$',
        ('pann', 'e_pm'): r'$\chi\chi \rightarrow e^-e^+$',
        ('gamma', 'gamma_gamma'): r'$\chi \rightarrow \gamma\gamma$',
        ('gamma', 'e_pm'): r'$\chi \rightarrow e^-e^+$',
    }

    xlims = {
        ('pann',  'gamma_gamma'): (1e-5, 5e3),
        ('pann',  'e_pm'):        (1e-3, 5e3),
        ('gamma', 'gamma_gamma'): (1e-5, 5e3),
        ('gamma', 'e_pm'):        (1e-3, 1e4),
    }

    for j, ch in enumerate(channels):
        ax = axes[0, j]
        xl_p = xlims[('pann', ch)]
        plot_panel(ax, 'pann', ch, xlim=xl_p, show_legend=False)
        ax.set_xlim(xl_p)
        ax.set_ylim(1e-29, 1e-25)
        ax.set_xlabel(r'$m_\chi[GeV]$', fontsize=FONT_SIZE)
        ax.set_ylabel(
            r'$\langle\sigma v\rangle/m_\chi\quad[cm^3s^{-1}GeV^{-1}]$',
            fontsize=FONT_SIZE,
        )
        ax.set_title(ch_titles[('pann', ch)], fontsize=TITLE_SIZE, pad=10)

        ax = axes[1, j]
        xl_g = xlims[('gamma', ch)]
        plot_panel(ax, 'gamma', ch, xlim=xl_g, show_legend=(j == 1))
        ax.set_xlim(xl_g)
        ax.set_ylim(1e-27, 1e-23)
        ax.set_xlabel(r'$m_\chi[GeV]$', fontsize=FONT_SIZE)
        ax.set_ylabel(r'$\Gamma_\chi[s^{-1}]$', fontsize=FONT_SIZE)
        ax.set_title(ch_titles[('gamma', ch)], fontsize=TITLE_SIZE, pad=10)

    fig.tight_layout(rect=[0, 0, 0.97, 0.96])

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig6_fisher.{fmt}')
        fig.savefig(path, dpi=dpi)
        print(f'Saved: {path}')

    for ch in channels:
        ref = PLANCK_REF.get(ch, {})
        suffix = '' if ch == 'gamma_gamma' else '_e_pm'
        raw_p = load_result(f'sig_pann_Planck{suffix}.npy')
        raw_g = load_result(f'sig_gamma_Planck{suffix}.npy')
        if raw_p is not None or raw_g is not None:
            print(f'\nPlanck calibration ({ch}):')
        if raw_p is not None:
            idx = np.argmin(np.abs(MASS_PANN - ref.get('ref_mass_GeV', 0.1)))
            ratio_p = ref.get('pann', 0) / raw_p[idx, 3] if raw_p[idx, 3] > 0 else 0
            print(f'  pann:  Fisher={raw_p[idx,3]:.2e} -> ref={ref.get("pann",0):.2e}  (x{ratio_p:.2f})')
        if raw_g is not None:
            idx = np.argmin(np.abs(MASS_GAMMA - ref.get('ref_mass_GeV', 0.1)))
            ratio_g = ref.get('gamma', 0) / raw_g[idx, 3] if raw_g[idx, 3] > 0 else 0
            print(f'  gamma: Fisher={raw_g[idx,3]:.2e} -> ref={ref.get("gamma",0):.2e}  (x{ratio_g:.2f})')


if __name__ == '__main__':
    plot_figure6()

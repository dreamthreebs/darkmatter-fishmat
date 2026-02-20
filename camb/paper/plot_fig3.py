#!/usr/bin/env python3
"""
Plot Figure 3 of arXiv:2304.07793:
  CMB power spectrum vs foreground residual vs instrumental noise.
  2x3 grid: rows = Ground / PICO, cols = TT / TE / EE.
  For TE: positive values solid, negative values dotted.
  TE has no instrumental noise (vanishes by construction).
  Noise plotted as black dashed; foreground as green (NILC V1).
"""

import sys
import os
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CAMB_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
FIGURE_DIR = os.path.join(SCRIPT_DIR, 'figures')
os.makedirs(FIGURE_DIR, exist_ok=True)

os.chdir(CAMB_DIR)
sys.path.insert(0, os.path.join(CAMB_DIR, 'calc'))

from consts import ells
from pdfunc import set_all_params, initial_totCL
from nls_fgres import ali_fg_res, pico_noise_level, pico_fg_res


def to_dl(cls_arr):
    ls = np.arange(len(cls_arr)) if cls_arr.ndim == 1 else np.arange(cls_arr.shape[0])
    fac = ls * (ls + 1) / (2 * np.pi)
    if cls_arr.ndim == 1:
        return fac * cls_arr
    return cls_arr * fac[:, np.newaxis]


# --- Analytical ALI noise (fitted to Ali_noise.dat) -----------------------
ALI_CHANNELS = [
    # (beam_arcmin, sigma_TT_deg, sigma_EE_deg)
    (19.0, 0.1352, 0.1904),   # 95 GHz
    (11.0, 0.2040, 0.2892),   # 150 GHz
]

def ali_noise_analytical(nells):
    ls = np.arange(nells, dtype=float)
    fac = ls * (ls + 1) / (2 * np.pi)
    fac[0] = 1e-30
    dl_out = np.zeros((nells, 2))
    for spec in range(2):
        inv_sum = np.zeros(nells)
        for beam_am, sig_TT, sig_EE in ALI_CHANNELS:
            sig = sig_TT if spec == 0 else sig_EE
            sig_rad = sig * 0.017453
            theta_rad = beam_am * 0.00029089
            Nl = sig_rad**2 * np.exp(ls * (ls + 1) * theta_rad**2 / (8 * np.log(2)))
            inv_sum += 1.0 / Nl
        Nl_comb = 1.0 / inv_sum
        dl_out[:, spec] = fac * Nl_comb
    return dl_out
# ---------------------------------------------------------------------------


def plot_signed_segments(ax, ls, dl, color, ls_pos='-', ls_neg=':', lw=1.2,
                         label=None):
    pos_done = False
    sign = np.sign(dl)
    i = 0
    n = len(dl)
    while i < n:
        if dl[i] == 0:
            i += 1
            continue
        s = sign[i]
        j = i
        while j < n and sign[j] == s:
            j += 1
        seg_l = ls[i:j]
        seg_d = np.abs(dl[i:j])
        if s > 0:
            lbl = label if not pos_done else None
            ax.loglog(seg_l, seg_d, color=color, ls=ls_pos, lw=lw, label=lbl)
            pos_done = True
        else:
            ax.loglog(seg_l, seg_d, color=color, ls=ls_neg, lw=lw)
        i = j


def plot_figure3():
    set_all_params()
    cmb_cls = initial_totCL()
    cmb_dl = to_dl(cmb_cls)

    ali_nls_dl = ali_noise_analytical(ells)
    ali_fg_dl = to_dl(ali_fg_res(ells))
    pico_nls_dl = to_dl(pico_noise_level(ells))
    pico_fg_dl = to_dl(pico_fg_res(ells))

    ls = np.arange(ells)

    # Axis ranges from the paper's actual figure
    xlims = {
        'Ground': (30, 1000),
        'PICO':   (10, 3000),
    }
    ylims = {
        ('Ground', 'TT'): (1e-3, 1e4),
        ('Ground', 'TE'): (1e-5, 1e2),
        ('Ground', 'EE'): (1e-2, 1e2),
        ('PICO',   'TT'): (1e-5, 1e4),
        ('PICO',   'TE'): (1e-6, 1e2),
        ('PICO',   'EE'): (1e-5, 1e2),
    }

    fg_color = 'tab:green'  # NILC V1 is green in the paper

    fig, axes = plt.subplots(2, 3, figsize=(11, 7))
    col_names = ['TT', 'TE', 'EE']
    spec_idx = [0, 2, 1]  # cls order: 0=TT, 1=EE, 2=TE

    experiments = [
        ('Ground', ali_fg_dl, ali_nls_dl),
        ('PICO',   pico_fg_dl, pico_nls_dl),
    ]

    for row, (exp_name, fg_dl, nls_dl) in enumerate(experiments):
        xl = xlims[exp_name]
        for j, (name, idx) in enumerate(zip(col_names, spec_idx)):
            ax = axes[row, j]
            mask = (ls >= xl[0]) & (ls <= xl[1])
            l_masked = ls[mask]

            cmb = cmb_dl[mask, idx]
            fg_raw = fg_dl[mask, idx]

            is_te = (name == 'TE')

            if is_te:
                plot_signed_segments(ax, l_masked, cmb, 'black',
                                     ls_pos='-', ls_neg=':', lw=1.6,
                                     label='CMB')
                valid_fg = fg_raw != 0
                if np.any(valid_fg):
                    plot_signed_segments(ax, l_masked[valid_fg], fg_raw[valid_fg],
                                         fg_color, ls_pos='-', ls_neg=':', lw=1.0,
                                         label='NILC, V1')
            else:
                nls_col = 0 if idx == 0 else 1
                noise = nls_dl[mask, nls_col]

                ax.loglog(l_masked, np.abs(cmb) + 1e-30,
                          'k-', lw=1.6, label='CMB')
                valid_fg = np.abs(fg_raw) > 0
                if np.any(valid_fg):
                    ax.loglog(l_masked[valid_fg], np.abs(fg_raw[valid_fg]),
                              color=fg_color, ls='-', lw=1.0, label='NILC, V1')
                valid_noise = noise > 0
                if np.any(valid_noise):
                    ax.loglog(l_masked[valid_noise], noise[valid_noise],
                              'k--', lw=1.2, label='Instrumental noise')

            ax.set_xlim(xl)
            ax.set_ylim(ylims[(exp_name, name)])
            ax.set_xlabel(r'$\ell$', fontsize=12)
            ax.grid(True, which='both', alpha=0.15)
            ax.tick_params(which='both', direction='in', top=True, right=True)

            ylabel_spec = {'TT': r'$D_\ell^{TT}$',
                           'TE': r'$D_\ell^{TE}$',
                           'EE': r'$D_\ell^{EE}$'}
            ax.set_ylabel(ylabel_spec[name] + r' [$\mu$K$^2$]', fontsize=11)

            # Legend in EE panels (right column), transparent background
            if j == 2:
                ax.legend(fontsize=7.5, loc='lower right', framealpha=0)

    fig.tight_layout()
    fig.subplots_adjust(hspace=0.4)

    for fmt, dpi in [('pdf', 300), ('png', 200)]:
        path = os.path.join(FIGURE_DIR, f'fig3_fg_compare.{fmt}')
        fig.savefig(path, dpi=dpi, bbox_inches='tight')
        print(f'Saved: {path}')


if __name__ == '__main__':
    plot_figure3()

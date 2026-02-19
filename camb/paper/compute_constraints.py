#!/usr/bin/env python3
"""
Compute Fisher matrix constraints on DM annihilation and decay parameters
for multiple experimental configurations.

Reproduces Figure 5 of arXiv:2304.07793 (gamma-gamma channel only).

Usage:
    conda activate test
    cd camb/paper
    python compute_constraints.py
"""

import sys
import os
import time
import numpy as np

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
CAMB_DIR = os.path.abspath(os.path.join(SCRIPT_DIR, '..'))
sys.path.insert(0, os.path.join(CAMB_DIR, 'calc'))
os.chdir(CAMB_DIR)

from consts import ells, params_num
from pdfunc import set_all_params, initial_totCL
from nls_fgres import (
    ali_noise_level, ali_fg_res,
    zero_noise, zero_fgres,
    fgres_dls2cls,
)
from sigfunc import (
    load_basic_6params_pd,
    get_combined_fisher_matrix,
)

RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')
os.makedirs(RESULTS_DIR, exist_ok=True)

# ============================================================
# Instrument specifications
# ============================================================

# PICO specs from arXiv:2304.07793 Table 1
# (freq_GHz, map_depth_P_uKarcmin, beam_FWHM_arcmin)
PICO_CHANNELS = [
    (21,   23.9,  38.4),
    (25,   18.4,  32.0),
    (30,   12.4,  28.3),
    (36,    7.9,  23.6),
    (43,    7.9,  22.2),
    (52,    5.7,  18.4),
    (62,    5.4,  12.8),
    (75,    4.2,  10.7),
    (90,    2.8,   9.5),
    (108,   2.3,   7.9),
    (129,   2.1,   7.4),
    (155,   1.8,   6.2),
    (186,   4.0,   4.3),
    (223,   4.5,   3.6),
    (268,   3.1,   3.2),
    (321,   4.2,   2.6),
    (385,   4.5,   2.5),
    (462,   9.1,   2.1),
    (555,  45.8,   1.5),
    (666, 177.0,   1.3),
    (799,1050.0,   1.1),
]

# Planck HFI+LFI specs (from Planck blue book / 2018 results)
# (freq_GHz, beam_FWHM_arcmin, noise_T_uKdeg, noise_P_uKdeg_or_None)
# 545/857 GHz have no polarization capability
PLANCK_CHANNELS = [
    (30,  32.29,  2.5,   3.5),
    (44,  27.94,  2.7,   4.0),
    (70,  13.08,  3.5,   5.0),
    (100,  9.66,  1.29,  1.96),
    (143,  7.22,  0.55,  1.17),
    (217,  4.90,  0.78,  1.75),
    (353,  4.92,  2.56,  7.31),
    (545,  4.67,  0.78,  None),
    (857,  4.22,  0.72,  None),
]


def _nls_single_channel(n_ells, beam_arcmin, depth_uKdeg):
    """N_l for a single frequency channel. depth in muK*deg, beam in arcmin."""
    arcmin2rad = np.pi / 180.0 / 60.0
    deg2rad = np.pi / 180.0
    ls = np.arange(n_ells)
    beam_rad = beam_arcmin * arcmin2rad
    sigma_rad = depth_uKdeg * deg2rad
    return sigma_rad**2 * np.exp(ls * (ls + 1) * beam_rad**2 / (8 * np.log(2)))


def pico_noise_analytical(n_ells):
    """Compute PICO combined noise N_l from 21 channels analytically.
    Returns array of shape (n_ells, 2) for [TT, EE] in C_l units."""
    arcmin2rad = np.pi / 180.0 / 60.0
    ls = np.arange(n_ells)

    inv_nls_TT = np.zeros(n_ells)
    inv_nls_EE = np.zeros(n_ells)

    for freq, depth_P, beam in PICO_CHANNELS:
        depth_T = depth_P / np.sqrt(2)
        beam_rad = beam * arcmin2rad

        exponent = ls * (ls + 1) * beam_rad**2 / (8 * np.log(2))
        beam_factor = np.exp(exponent)

        nls_T = (depth_T * arcmin2rad)**2 * beam_factor
        nls_E = (depth_P * arcmin2rad)**2 * beam_factor

        with np.errstate(divide='ignore'):
            inv_nls_TT += 1.0 / nls_T
            inv_nls_EE += 1.0 / nls_E

    with np.errstate(divide='ignore', invalid='ignore'):
        nls_TT = np.where(inv_nls_TT > 0, 1.0 / inv_nls_TT, 0.0)
        nls_EE = np.where(inv_nls_EE > 0, 1.0 / inv_nls_EE, 0.0)

    return np.stack((nls_TT, nls_EE), axis=1)


def planck_noise_analytical(n_ells):
    """Compute Planck combined noise N_l analytically.
    545/857 GHz channels contribute only to TT (no polarization).
    Returns array of shape (n_ells, 2) for [TT, EE] in C_l units."""
    inv_nls_TT = np.zeros(n_ells)
    inv_nls_EE = np.zeros(n_ells)

    for freq, beam, noise_T, noise_P in PLANCK_CHANNELS:
        nls_T = _nls_single_channel(n_ells, beam, noise_T)
        with np.errstate(divide='ignore'):
            inv_nls_TT += 1.0 / nls_T

        if noise_P is not None:
            nls_E = _nls_single_channel(n_ells, beam, noise_P)
            with np.errstate(divide='ignore'):
                inv_nls_EE += 1.0 / nls_E

    with np.errstate(divide='ignore', invalid='ignore'):
        nls_TT = np.where(inv_nls_TT > 0, 1.0 / inv_nls_TT, 0.0)
        nls_EE = np.where(inv_nls_EE > 0, 1.0 / inv_nls_EE, 0.0)

    return np.stack((nls_TT, nls_EE), axis=1)


def pico_fgres_from_file(n_ells):
    """Load PICO foreground residual from ILC simulation files.
    Data covers l up to ~1024; beyond that we pad with zeros."""
    TT = np.load('./input/pico/fg_TT_results.npy')
    EE = np.load('./input/pico/fg_EE_results.npy')
    TE = np.load('./input/pico/fg_TE_results.npy')
    half = np.stack((TT, EE, TE), axis=1)
    pad = np.zeros((n_ells - len(half), 3))
    fgres_dls = np.insert(half, len(half), pad, axis=0)
    return fgres_dls2cls(fgres_dls, n_ells, 1, len(TT))


# ============================================================
# Fisher matrix computation for a batch of DM masses
# ============================================================

def compute_2sigma_vs_mass(pd_dm_diff_mass, cls, nls, fgres, lmin, lmax, fsky):
    """Compute 2-sigma Fisher constraints for each DM mass point.

    Args:
        pd_dm_diff_mass: shape (ells, 3, n_mass) — dC_l/d(DM param) at each mass
        cls: shape (ells, 3) — fiducial C_l [TT, EE, TE]
        nls: shape (ells, 2) — noise [TT, EE]
        fgres: shape (ells, 3) — foreground residual [TT, EE, TE]
        lmin, lmax: multipole range
        fsky: sky fraction

    Returns:
        sigma_2: shape (n_mass, 4) — 2-sigma bounds [TT, EE, TE, combined]
    """
    (thetastarmc_pd, ombh2_pd, omch2_pd,
     As_pd, ns_pd, optical_depth_pd) = load_basic_6params_pd()

    pd = np.zeros((ells, 3, params_num))
    pd[:, :, 0] = ombh2_pd
    pd[:, :, 1] = omch2_pd
    pd[:, :, 2] = thetastarmc_pd
    pd[:, :, 3] = optical_depth_pd
    pd[:, :, 4] = ns_pd
    pd[:, :, 5] = As_pd

    n_mass = pd_dm_diff_mass.shape[2]
    sigma_2 = np.zeros((n_mass, 4))

    for i in range(n_mass):
        pd[:, :, 6] = pd_dm_diff_mass[:, :, i]
        sig_combined = get_combined_fisher_matrix(
            pd, cls, nls, fgres, lmin, lmax, fsky
        )
        sigma_2[i, 3] = float(2 * sig_combined[params_num - 1])

        if (i + 1) % 10 == 0 or i == 0:
            print(f'  mass point {i+1}/{n_mass}, '
                  f'2sigma(combined) = {sigma_2[i,3]:.3e}')

    return sigma_2


# ============================================================
# Experimental configurations matching the paper
# ============================================================

CONFIGS = {
    'CVL': {
        'lmin': 10, 'lmax': 3000, 'fsky': 1.0,
        'noise_func': lambda: zero_noise(ells),
        'fgres_func': lambda: zero_fgres(ells),
        'description': 'Cosmic Variance Limited, l=[10,3000], fsky=1',
    },
    'ALI_fgres': {
        'lmin': 30, 'lmax': 620, 'fsky': 0.4,
        'noise_func': lambda: ali_noise_level(ells),
        'fgres_func': lambda: ali_fg_res(ells),
        'description': 'Ground/ALI, l=[30,620], fsky=0.4, with foreground',
    },
    'ALI_nofgres': {
        'lmin': 30, 'lmax': 620, 'fsky': 0.4,
        'noise_func': lambda: ali_noise_level(ells),
        'fgres_func': lambda: zero_fgres(ells),
        'description': 'Ground/ALI, l=[30,620], fsky=0.4, noise only',
    },
    'PICO_fgres': {
        'lmin': 10, 'lmax': 3000, 'fsky': 0.77,
        'noise_func': lambda: pico_noise_analytical(ells),
        'fgres_func': lambda: pico_fgres_from_file(ells),
        'description': 'PICO, l=[10,3000], fsky=0.77, with foreground',
    },
    'PICO_nofgres': {
        'lmin': 10, 'lmax': 3000, 'fsky': 0.77,
        'noise_func': lambda: pico_noise_analytical(ells),
        'fgres_func': lambda: zero_fgres(ells),
        'description': 'PICO, l=[10,3000], fsky=0.77, noise only',
    },
    'Planck': {
        'lmin': 30, 'lmax': 2500, 'fsky': 0.57,
        'noise_func': lambda: planck_noise_analytical(ells),
        'fgres_func': lambda: zero_fgres(ells),
        'description': 'Planck (Fisher, ideal noise), l=[30,2500], fsky=0.57',
    },
}


def main():
    t0 = time.time()

    print('=' * 60)
    print('Fisher Matrix Constraints — arXiv:2304.07793 Figure 5')
    print('Channel: DM → γγ')
    print('=' * 60)

    # Load pre-computed partial derivatives (50 mass points)
    print('\nLoading pre-computed partial derivatives...')
    pd_pann = np.load('./data/pd_pann_50_data/pd_pann_50_diff_mass.npy')
    pd_gamma = np.load('./data/pd_gamma_50_data/pd_gamma_50_diff_mass.npy')
    print(f'  pd_pann shape: {pd_pann.shape}')
    print(f'  pd_gamma shape: {pd_gamma.shape}')

    # Compute fiducial C_l (runs CAMB once)
    print('\nComputing fiducial C_l (running CAMB)...')
    set_all_params()
    cls = initial_totCL()
    print(f'  cls shape: {cls.shape}')

    # Loop over configurations
    for config_name, cfg in CONFIGS.items():
        print(f'\n{"="*60}')
        print(f'Config: {config_name}')
        print(f'  {cfg["description"]}')
        print(f'{"="*60}')

        nls = cfg['noise_func']()
        fgres = cfg['fgres_func']()

        # --- Annihilation (p_ann) ---
        print(f'\n  Computing p_ann constraints...')
        t1 = time.time()
        sig_pann = compute_2sigma_vs_mass(
            pd_pann, cls, nls, fgres,
            cfg['lmin'], cfg['lmax'], cfg['fsky']
        )
        fname = os.path.join(RESULTS_DIR, f'sig_pann_{config_name}.npy')
        np.save(fname, sig_pann)
        print(f'  Saved: {fname} ({time.time()-t1:.1f}s)')

        # --- Decay (Gamma) ---
        print(f'\n  Computing Gamma constraints...')
        t1 = time.time()
        sig_gamma = compute_2sigma_vs_mass(
            pd_gamma, cls, nls, fgres,
            cfg['lmin'], cfg['lmax'], cfg['fsky']
        )
        fname = os.path.join(RESULTS_DIR, f'sig_gamma_{config_name}.npy')
        np.save(fname, sig_gamma)
        print(f'  Saved: {fname} ({time.time()-t1:.1f}s)')

    print(f'\nAll done! Total time: {time.time()-t0:.1f}s')
    print(f'Results saved to: {RESULTS_DIR}/')


if __name__ == '__main__':
    main()

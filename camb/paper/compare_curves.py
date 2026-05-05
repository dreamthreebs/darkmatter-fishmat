#!/usr/bin/env python3
"""
专门回答两个问题：
  Fig 4: 每个子图里，黑线(TT) 比 红线(EE) 弱多少 ；蓝线(TT+TE+EE) 比 红线(EE) 低多少。
  Fig 5: 每个子图里，加了前景(实线) 比 没加前景(虚线) 的约束退化多少。

定义 (合理且符合 log-log 图惯例)：
  曲线 y(m) 是 2σ 上限，y 越小代表约束越强。
  比值     r(m) = y_A(m) / y_B(m)
  "A 比 B 弱 X 倍" ⇔ r > 1，报告 X = r；"A 比 B 低/好 Y 倍" ⇔ r < 1，报告 Y = 1/r。
  因为上限在质量轴上跨多个数量级，我们用 **对数平均(几何平均)** 作为全局代表值：
      G(r) = exp( mean(ln r) )
  同时报告质量区间内的最小/最大比值，以及五个质量点(0.01,0.1,1,10,100 GeV)的逐点值。

Usage:
    cd camb/paper
    python compare_curves.py
"""

import os
import numpy as np
from scipy.interpolate import CubicSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')

MASS_PANN  = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)


def load(name):
    path = os.path.join(RESULTS_DIR, name)
    return np.load(path) if os.path.exists(path) else None


def interp_at(mass_grid, y, m):
    valid = y > 0
    if np.sum(valid) < 4:
        return np.nan
    cs = CubicSpline(np.log10(mass_grid[valid]), np.log10(y[valid]))
    return 10 ** cs(np.log10(m))


def geo_mean(r):
    r = r[np.isfinite(r) & (r > 0)]
    if r.size == 0:
        return np.nan
    return np.exp(np.mean(np.log(r)))


def summarize_ratio(label_A_over_B, y_A, y_B, mass_grid, mass_points):
    """Return dict describing ratio y_A / y_B."""
    both_valid = (y_A > 0) & (y_B > 0) & np.isfinite(y_A) & np.isfinite(y_B)
    r = y_A[both_valid] / y_B[both_valid]
    g = geo_mean(r)
    rmin, rmax = (np.min(r), np.max(r)) if r.size else (np.nan, np.nan)
    print(f'  {label_A_over_B}:')
    print(f'    几何平均 = {g:.2f} x    '
          f'范围 [{rmin:.2f} x , {rmax:.2f} x ]   (N={r.size} 质量点)')
    row = '    逐点: '
    for m in mass_points:
        yA = interp_at(mass_grid, y_A, m)
        yB = interp_at(mass_grid, y_B, m)
        row += f'm={m:<5g} -> {yA/yB:5.2f} x ;   '
    print(row.rstrip(';   '))
    return g, rmin, rmax


def sec(title):
    print('\n' + '=' * 72)
    print('  ' + title)
    print('=' * 72)


def fig4():
    """Fig 4: ALI (Ground), γγ, TT=col0 (black), EE=col1 (red), TT+TE+EE=col3 (blue)."""
    sec('Fig 4 — Ground (ALI), γγ 通道，比较 TT / EE / TT+TE+EE')
    print('  列: 0=TT(黑)  1=EE(红)  3=TT+TE+EE(蓝)')
    print('  数值越小代表约束越强。比值 = (y_A / y_B)。')

    mass_points_pann  = [0.01, 0.1, 1, 10, 100]
    mass_points_gamma = [0.01, 0.1, 1, 10, 100]

    for pt, label, mass_grid, mps in [
        ('pann',  r'湮灭 ⟨σv⟩/m_χ', MASS_PANN,  mass_points_pann),
        ('gamma', r'衰变 Γ_χ',     MASS_GAMMA, mass_points_gamma),
    ]:
        d = load(f'sig_{pt}_ALI_fgres.npy')
        if d is None:
            print(f'\n[缺] sig_{pt}_ALI_fgres.npy')
            continue

        tt, ee, comb = d[:, 0], d[:, 1], d[:, 3]

        print(f'\n-- {label} --')
        summarize_ratio('TT / EE    (黑/红, >1 表示黑线比红线弱)',
                        tt, ee, mass_grid, mps)
        summarize_ratio('TT+TE+EE / EE    (蓝/红, <1 表示蓝线比红线低/更强)',
                        comb, ee, mass_grid, mps)


def fig5():
    """Fig 5: fgres vs nofgres for Ground (ALI, red) and PICO (blue)."""
    sec('Fig 5 — 每个子图里 加前景(实线) / 不加前景(虚线) 的约束退化')
    print('  比值 = y_fg / y_nofg;  >1 表示加了前景后上限变大(约束变弱)。')
    print('  全局统计: 50 个质量点上 log 域的几何平均 + 极值。')

    experiments = [
        ('ALI',  'Ground (ALI, red)'),
        ('PICO', 'PICO (blue)'),
    ]
    params = [
        ('pann',  r'湮灭 ⟨σv⟩/m_χ', MASS_PANN),
        ('gamma', r'衰变 Γ_χ',     MASS_GAMMA),
    ]
    mass_points = [0.01, 0.1, 1, 10, 100]

    for exp_name, exp_label in experiments:
        print(f'\n### {exp_label}')
        for pt, plabel, mass_grid in params:
            d_fg   = load(f'sig_{pt}_{exp_name}_fgres.npy')
            d_nofg = load(f'sig_{pt}_{exp_name}_nofgres.npy')
            if d_fg is None or d_nofg is None:
                print(f'  [缺] {exp_name} {pt}')
                continue
            print(f'\n-- {plabel} --')
            summarize_ratio('fgres / nofgres  (实线/虚线)',
                            d_fg[:, 3], d_nofg[:, 3], mass_grid, mass_points)

            r = d_fg[:, 3] / d_nofg[:, 3]
            valid = (d_fg[:, 3] > 0) & (d_nofg[:, 3] > 0)
            pct = (r[valid] - 1) * 100
            print(f'    退化百分比:  几何平均 {geo_mean(r[valid]) * 100 - 100:+.1f}%   '
                  f'算术平均 {np.mean(pct):+.1f}%   '
                  f'范围 [{np.min(pct):+.1f}%, {np.max(pct):+.1f}%]')


def fig6():
    """Fig 6: Planck (calibrated) vs Ground (ALI) vs PICO vs CVL,
       improvements at m=0.1 GeV for both γγ and e+e- channels.
    """
    sec('Fig 6 — m_χ = 0.1 GeV 处各实验对 Planck 的改善倍数')
    print('  Planck 曲线: Fisher 预测 (sig_*_Planck*.npy 第 4 列), 在 m=0.1 GeV 处')
    print('  按文献值校准——γγ pann -> 6.0e-28, γγ gamma -> 2.5e-25,')
    print('  e+e- pann -> 4.5e-28, e+e- gamma -> 5.0e-26 (cm^3 s^-1 GeV^-1 / s^-1)')
    print('  曲线形状沿用 Fisher 预测、整体仅做一个常数 rescale。')
    print('  "A 比 B 改善 k 倍" 定义为 y_B / y_A = k (>1，A 给出更紧的上限)。')

    PLANCK_REF = {
        'gamma_gamma': {'pann': 6.0e-28, 'gamma': 2.5e-25, 'ref_mass_GeV': 0.1},
        'e_pm':        {'pann': 4.5e-28, 'gamma': 5.0e-26, 'ref_mass_GeV': 0.1},
    }

    M_REF = 0.1
    for ch_label, suffix, ch_show in [
        ('gamma_gamma', '', 'γγ'),
        ('e_pm',        '_e_pm', 'e+e-'),
    ]:
        ref = PLANCK_REF[ch_label]
        print(f'\n### 通道 {ch_show}')
        for pt, p_label, mass_grid, unit in [
            ('pann',  '湮灭 ⟨σv⟩/m_χ', MASS_PANN,  'cm^3 s^-1 GeV^-1'),
            ('gamma', '衰变 Γ_χ',     MASS_GAMMA, 's^-1'),
        ]:
            raw = load(f'sig_{pt}_Planck{suffix}.npy')
            if raw is None:
                continue
            idx = np.argmin(np.abs(mass_grid - ref['ref_mass_GeV']))
            cal_factor = ref[pt] / raw[idx, 3]
            planck_y = raw[:, 3] * cal_factor

            ali  = load(f'sig_{pt}_ALI_fgres{suffix}.npy')
            pico = load(f'sig_{pt}_PICO_fgres{suffix}.npy')
            cvl  = load(f'sig_{pt}_CVL{suffix}.npy')

            p   = interp_at(mass_grid, planck_y, M_REF)
            ali_v  = interp_at(mass_grid, ali[:, 3],  M_REF) if ali  is not None else np.nan
            pico_v = interp_at(mass_grid, pico[:, 3], M_REF) if pico is not None else np.nan
            cvl_v  = interp_at(mass_grid, cvl[:, 3],  M_REF) if cvl  is not None else np.nan

            print(f'\n-- {p_label}  [{unit}], m_χ=0.1 GeV --')
            print(f'    Planck (校准)        = {p:.2e}')
            print(f'    Ground (ALI+fg)     = {ali_v:.2e}      Planck/Ground = {p/ali_v:.2f}x')
            print(f'    PICO (PICO+fg)      = {pico_v:.2e}      Planck/PICO   = {p/pico_v:.2f}x'
                  f'   Ground/PICO = {ali_v/pico_v:.2f}x')
            print(f'    CVL                 = {cvl_v:.2e}      Planck/CVL    = {p/cvl_v:.2f}x'
                  f'   PICO/CVL    = {pico_v/cvl_v:.2f}x')

            print(f'    论文写法:')
            print(f'      "对于{p_label.split()[0]}通道, 以 m_χ ~ 0.1 GeV、{ch_show} 为例, '
                  f'地面观测对{p_label.split()[0]}参数的 95% 置信上限达到 '
                  f'{p_label.split()[1]} < {ali_v:.1e} {unit},')
            print(f'       相对于 Planck 改善约 {p/ali_v:.1f} 倍, '
                  f'PICO 在此基础上可再进一步改善约 {ali_v/pico_v:.1f} 倍."')


def main():
    print('定义:')
    print('  - 2σ 上限曲线 y(m),  y 越小 = 约束越强。')
    print('  - 比值 r = y_A / y_B.  "A 比 B 弱 k 倍" 意为 r = k (>1).')
    print('  - 全局指标用 log 域几何平均（因为 y 在 m 上跨多个数量级）。')
    fig4()
    fig5()
    fig6()


if __name__ == '__main__':
    main()

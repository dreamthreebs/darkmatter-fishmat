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


def main():
    print('定义:')
    print('  - 2σ 上限曲线 y(m),  y 越小 = 约束越强。')
    print('  - 比值 r = y_A / y_B.  "A 比 B 弱 k 倍" 意为 r = k (>1).')
    print('  - 全局指标用 log 域几何平均（因为 y 在 m 上跨多个数量级）。')
    fig4()
    fig5()


if __name__ == '__main__':
    main()

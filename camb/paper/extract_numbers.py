#!/usr/bin/env python3
"""
从 results/*.npy 自动提取论文中需要引用的所有数值。
换 fg residual 后重跑此脚本即可拿到所有要写进正文的数字。

Usage:
    conda activate test
    cd camb/paper
    python extract_numbers.py
"""

import os
import numpy as np
from scipy.interpolate import CubicSpline

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
RESULTS_DIR = os.path.join(SCRIPT_DIR, 'results')

MASS_PANN  = np.geomspace(1e-5, 5e3, 50)
MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)

PLANCK_REF = {
    'gamma_gamma': {'pann': 6.0e-28, 'gamma': 2.5e-25, 'ref_mass_GeV': 0.1},
    'e_pm':        {'pann': 4.5e-28, 'gamma': 5.0e-26, 'ref_mass_GeV': 0.1},
}

THERMAL_SV = 3e-26  # cm^3/s


def load(name):
    path = os.path.join(RESULTS_DIR, name)
    if os.path.exists(path):
        return np.load(path)
    return None


def interp_at_mass(mass_grid, data_col, m_target):
    """log-log cubic interpolation, return value at m_target."""
    valid = data_col > 0
    if np.sum(valid) < 4:
        return np.nan
    cs = CubicSpline(np.log10(mass_grid[valid]), np.log10(data_col[valid]))
    return 10**cs(np.log10(m_target))


def calibrate_planck(fisher_data, mass_grid, ref_value, ref_mass):
    idx = np.argmin(np.abs(mass_grid - ref_mass))
    ratio = ref_value / fisher_data[idx, 3] if fisher_data[idx, 3] > 0 else 1.0
    out = fisher_data.copy()
    out[:, 3] = fisher_data[:, 3] * ratio
    return out


def find_thermal_wimp_mass(mass_grid, pann_combined):
    """Find mass where 2σ(p_ann) crosses thermal relic line 3e-26/m_χ.
    Returns m_χ* below which thermal WIMPs are excluded (constraint < thermal)."""
    thermal = THERMAL_SV / mass_grid
    ratio = pann_combined / thermal
    log_ratio = np.log10(ratio)
    for i in range(len(log_ratio) - 1):
        if log_ratio[i] <= 0 and log_ratio[i+1] > 0:
            frac = -log_ratio[i] / (log_ratio[i+1] - log_ratio[i])
            m_cross = 10**(np.log10(mass_grid[i]) + frac *
                           (np.log10(mass_grid[i+1]) - np.log10(mass_grid[i])))
            return m_cross
    return np.nan


def sec(title):
    print(f'\n{"="*70}')
    print(f'  {title}')
    print(f'{"="*70}')


def subsec(title):
    print(f'\n--- {title} ---')


def main():
    print('论文数值提取 — 所有需要写入正文的数字')
    print('='*70)

    # ================================================================
    # 1) Fig 4: TT vs EE vs TT+TE+EE 比值 (Ground, γγ)
    # ================================================================
    sec('1. Fig 4: TT vs EE vs Combined 比值 (Ground/ALI, γγ)')

    for pt, label, mass_grid in [
        ('pann', '湮灭 ⟨σv⟩/m_χ', MASS_PANN),
        ('gamma', '衰变 Γ_χ', MASS_GAMMA),
    ]:
        d = load(f'sig_{pt}_ALI_fgres.npy')
        if d is None:
            print(f'  [缺] sig_{pt}_ALI_fgres.npy')
            continue

        subsec(f'{label}')
        tt = d[:, 0]
        ee = d[:, 1]
        combined = d[:, 3]

        ratio_tt_ee = tt / ee
        ratio_ee_comb = ee / combined

        avg_tt_ee = np.nanmean(ratio_tt_ee[ratio_tt_ee > 0])
        med_tt_ee = np.nanmedian(ratio_tt_ee[ratio_tt_ee > 0])
        avg_ee_comb = np.nanmean(ratio_ee_comb[ratio_ee_comb > 0])
        med_ee_comb = np.nanmedian(ratio_ee_comb[ratio_ee_comb > 0])

        print(f'  TT/EE  (EE比TT严多少倍): 平均={avg_tt_ee:.1f}x, 中位数={med_tt_ee:.1f}x')
        print(f'  EE/Combined (Combined比EE好多少倍): 平均={avg_ee_comb:.1f}x, 中位数={med_ee_comb:.1f}x')

        for m in [0.1, 1, 10, 100, 1000]:
            tt_m = interp_at_mass(mass_grid, tt, m)
            ee_m = interp_at_mass(mass_grid, ee, m)
            cb_m = interp_at_mass(mass_grid, combined, m)
            if not np.isnan(cb_m):
                print(f'  m={m:>6g} GeV:  TT={tt_m:.2e}  EE={ee_m:.2e}  Combined={cb_m:.2e}  '
                      f'TT/EE={tt_m/ee_m:.1f}x  EE/Comb={ee_m/cb_m:.1f}x')

    # ================================================================
    # 2) Fig 5: 前景残差 degradation (noise+fg vs noise-only)
    # ================================================================
    sec('2. Fig 5: 前景残差导致的约束退化 (noise+fg vs noise-only)')

    for channel, suffix in [('γγ', ''), ('e⁺e⁻', '_e_pm')]:
        subsec(f'通道: {channel}')
        for exp_name, exp_label in [('ALI', 'Ground'), ('PICO', 'PICO')]:
            for pt, pt_label, mass_grid in [
                ('pann', '湮灭', MASS_PANN),
                ('gamma', '衰变', MASS_GAMMA),
            ]:
                d_fg   = load(f'sig_{pt}_{exp_name}_fgres{suffix}.npy')
                d_nofg = load(f'sig_{pt}_{exp_name}_nofgres{suffix}.npy')
                if d_fg is None or d_nofg is None:
                    continue

                ratio = d_fg[:, 3] / d_nofg[:, 3]
                valid = (d_fg[:, 3] > 0) & (d_nofg[:, 3] > 0)
                if np.sum(valid) == 0:
                    continue

                degradation_pct = (ratio[valid] - 1) * 100
                avg_deg = np.nanmean(degradation_pct)
                min_deg = np.nanmin(degradation_pct)
                max_deg = np.nanmax(degradation_pct)

                print(f'  {exp_label} {pt_label}: '
                      f'平均退化 {avg_deg:+.1f}%, 范围 [{min_deg:+.1f}%, {max_deg:+.1f}%]')

    # ================================================================
    # 3) Fig 6: 特定质量点的约束值 (m_χ = 0.1 GeV 等)
    # ================================================================
    sec('3. Fig 6: m_χ = 0.1 GeV 处各实验的约束 (95% C.L.)')

    for channel, suffix, ch_label in [('gamma_gamma', '', 'γγ'), ('e_pm', '_e_pm', 'e⁺e⁻')]:
        subsec(f'通道: {ch_label}')
        ref = PLANCK_REF[channel]
        m_ref = ref['ref_mass_GeV']

        for pt, pt_label, mass_grid, unit in [
            ('pann', '⟨σv⟩/m_χ', MASS_PANN, 'cm³ s⁻¹ GeV⁻¹'),
            ('gamma', 'Γ_χ', MASS_GAMMA, 's⁻¹'),
        ]:
            print(f'\n  {pt_label} at m_χ={m_ref} GeV [{unit}]:')

            planck_raw = load(f'sig_{pt}_Planck{suffix}.npy')
            planck_val = np.nan
            if planck_raw is not None:
                planck_cal = calibrate_planck(planck_raw, mass_grid,
                                              ref[pt if pt in ref else 'gamma'],
                                              ref['ref_mass_GeV'])
                planck_val = interp_at_mass(mass_grid, planck_cal[:, 3], m_ref)
                print(f'    Planck (校准后)   = {planck_val:.2e}')

            for exp_name, exp_label in [
                ('ALI_fgres', 'Ground (noise+fg)'),
                ('ALI_nofgres', 'Ground (noise only)'),
                ('PICO_fgres', 'PICO (noise+fg)'),
                ('PICO_nofgres', 'PICO (noise only)'),
                ('CVL', 'CVL'),
            ]:
                d = load(f'sig_{pt}_{exp_name}{suffix}.npy')
                if d is None:
                    continue
                val = interp_at_mass(mass_grid, d[:, 3], m_ref)
                ratio_str = ''
                if not np.isnan(planck_val) and planck_val > 0 and val > 0:
                    ratio_str = f'  (Planck/this = {planck_val/val:.1f}x)'
                print(f'    {exp_label:25s} = {val:.2e}{ratio_str}')

    # ================================================================
    # 4) Planck 校准系数
    # ================================================================
    sec('4. Planck 校准系数')

    for channel, suffix, ch_label in [('gamma_gamma', '', 'γγ'), ('e_pm', '_e_pm', 'e⁺e⁻')]:
        ref = PLANCK_REF[channel]
        print(f'\n  通道 {ch_label} (参考质量 = {ref["ref_mass_GeV"]} GeV):')
        for pt, pt_label, mass_grid in [
            ('pann', 'pann', MASS_PANN),
            ('gamma', 'Gamma', MASS_GAMMA),
        ]:
            raw = load(f'sig_{pt}_Planck{suffix}.npy')
            if raw is None:
                print(f'    {pt_label}: [缺文件]')
                continue
            idx = np.argmin(np.abs(mass_grid - ref['ref_mass_GeV']))
            fisher_val = raw[idx, 3]
            ref_val = ref[pt if pt in ref else 'gamma']
            ratio = ref_val / fisher_val if fisher_val > 0 else 0
            print(f'    {pt_label}: Fisher={fisher_val:.2e} -> ref={ref_val:.2e}  (x{ratio:.2f})')

    # ================================================================
    # 5) 各实验之间的倍数比较
    # ================================================================
    sec('5. 实验间约束的倍数比较 (m_χ = 0.1 GeV)')

    for channel, suffix, ch_label in [('gamma_gamma', '', 'γγ'), ('e_pm', '_e_pm', 'e⁺e⁻')]:
        subsec(f'通道: {ch_label}')
        ref = PLANCK_REF[channel]
        m_ref = ref['ref_mass_GeV']

        for pt, pt_label, mass_grid in [
            ('pann', '湮灭', MASS_PANN),
            ('gamma', '衰变', MASS_GAMMA),
        ]:
            planck_raw = load(f'sig_{pt}_Planck{suffix}.npy')
            if planck_raw is None:
                continue
            planck_cal = calibrate_planck(planck_raw, mass_grid,
                                          ref[pt if pt in ref else 'gamma'],
                                          ref['ref_mass_GeV'])
            p_val = interp_at_mass(mass_grid, planck_cal[:, 3], m_ref)

            exps = {}
            for name, label in [('ALI_fgres', 'Ground'), ('PICO_fgres', 'PICO'), ('CVL', 'CVL')]:
                d = load(f'sig_{pt}_{name}{suffix}.npy')
                if d is not None:
                    exps[label] = interp_at_mass(mass_grid, d[:, 3], m_ref)

            print(f'\n  {pt_label}:')
            print(f'    Planck (校准后) = {p_val:.2e}')
            for label, val in exps.items():
                ratio = p_val / val if val > 0 else np.nan
                print(f'    {label:10s} = {val:.2e}   Planck比{label}宽 {ratio:.1f}x')

            if 'Ground' in exps and 'PICO' in exps:
                print(f'    PICO比Ground好 {exps["Ground"]/exps["PICO"]:.1f}x')
            if 'PICO' in exps and 'CVL' in exps:
                print(f'    CVL比PICO好 {exps["PICO"]/exps["CVL"]:.1f}x')

    # ================================================================
    # 6) Thermal WIMP 排除质量 (热遗迹暗物质)
    # ================================================================
    sec('6. Thermal WIMP 排除质量阈值 (⟨σv⟩ = 3×10⁻²⁶ cm³/s)')
    print('  thermal relic line: p_ann = 3e-26/m_χ')
    print('  交叉点 m_χ*: 约束曲线与热遗迹线交点，m < m_χ* 被排除')

    for channel, suffix, ch_label in [('gamma_gamma', '', 'γγ'), ('e_pm', '_e_pm', 'e⁺e⁻')]:
        subsec(f'通道: {ch_label}')
        ref = PLANCK_REF[channel]

        planck_raw = load(f'sig_pann_Planck{suffix}.npy')
        if planck_raw is not None:
            planck_cal = calibrate_planck(planck_raw, MASS_PANN,
                                          ref['pann'], ref['ref_mass_GeV'])
            m_cross = find_thermal_wimp_mass(MASS_PANN, planck_cal[:, 3])
            if np.isnan(m_cross):
                print(f'  {"Planck (校准后)":25s}: 未找到交叉点')
            else:
                print(f'  {"Planck (校准后)":25s}: m_χ < {m_cross:.0f} GeV 被排除')

        for exp_name, exp_label in [
            ('ALI_fgres', 'Ground (noise+fg)'),
            ('ALI_nofgres', 'Ground (noise only)'),
            ('PICO_fgres', 'PICO (noise+fg)'),
            ('PICO_nofgres', 'PICO (noise only)'),
            ('CVL', 'CVL'),
        ]:
            d = load(f'sig_pann_{exp_name}{suffix}.npy')
            if d is None:
                continue
            m_cross = find_thermal_wimp_mass(MASS_PANN, d[:, 3])
            if np.isnan(m_cross):
                print(f'  {exp_label:25s}: 未找到交叉点（约束不够强）')
            else:
                print(f'  {exp_label:25s}: m_χ < {m_cross:.0f} GeV 被排除')

    # ================================================================
    # 7) 前景退化的全局平均 (论文结论部分)
    # ================================================================
    sec('7. 前景退化全局汇总 (跨质量点、跨通道)')

    all_degradations = []
    for channel, suffix in [('gamma_gamma', ''), ('e_pm', '_e_pm')]:
        for exp_name in ['ALI', 'PICO']:
            for pt, mass_grid in [('pann', MASS_PANN), ('gamma', MASS_GAMMA)]:
                d_fg   = load(f'sig_{pt}_{exp_name}_fgres{suffix}.npy')
                d_nofg = load(f'sig_{pt}_{exp_name}_nofgres{suffix}.npy')
                if d_fg is None or d_nofg is None:
                    continue
                ratio = d_fg[:, 3] / d_nofg[:, 3]
                valid = (d_fg[:, 3] > 0) & (d_nofg[:, 3] > 0)
                if np.sum(valid) > 0:
                    degs = (ratio[valid] - 1) * 100
                    all_degradations.extend(degs.tolist())

    if all_degradations:
        arr = np.array(all_degradations)
        print(f'  所有实验+通道+参数的前景退化:')
        print(f'    范围: {np.min(arr):+.1f}% ~ {np.max(arr):+.1f}%')
        print(f'    平均: {np.mean(arr):+.1f}%')
        print(f'  论文写法: "CMB前景残余使约束减弱 {np.min(arr):.0f}%–{np.max(arr):.0f}%"')

    # ================================================================
    # 8) 更多质量点的约束值 (可做完整表格)
    # ================================================================
    sec('8. 完整表格: 各实验在多个质量点的 Combined 2σ 上限')

    mass_points = [0.01, 0.1, 1, 10, 100, 1000]

    for channel, suffix, ch_label in [('gamma_gamma', '', 'γγ'), ('e_pm', '_e_pm', 'e⁺e⁻')]:
        for pt, pt_label, mass_grid, unit in [
            ('pann', '⟨σv⟩/m_χ', MASS_PANN, 'cm³ s⁻¹ GeV⁻¹'),
            ('gamma', 'Γ_χ', MASS_GAMMA, 's⁻¹'),
        ]:
            subsec(f'{ch_label}, {pt_label} [{unit}]')

            header = f'{"m_χ [GeV]":>12s}'
            exp_list = ['ALI_fgres', 'ALI_nofgres', 'PICO_fgres', 'PICO_nofgres', 'CVL']
            exp_labels = ['Ground+fg', 'Ground', 'PICO+fg', 'PICO', 'CVL']

            exp_data = {}
            for name in exp_list:
                d = load(f'sig_{pt}_{name}{suffix}.npy')
                exp_data[name] = d

            for lbl in exp_labels:
                header += f'  {lbl:>12s}'
            print(header)
            print('-' * len(header))

            for m in mass_points:
                if m < mass_grid[0] or m > mass_grid[-1]:
                    continue
                row = f'{m:>12g}'
                for name in exp_list:
                    d = exp_data[name]
                    if d is not None:
                        val = interp_at_mass(mass_grid, d[:, 3], m)
                        row += f'  {val:>12.2e}'
                    else:
                        row += f'  {"—":>12s}'
                print(row)

    print('\n' + '='*70)
    print('  完成。以上所有数值可直接用于论文正文。')
    print('  换 fg 后重跑: python compute_constraints.py && python extract_numbers.py')
    print('='*70)


if __name__ == '__main__':
    main()

# Fig 4 / Fig 5 曲线比值汇总

本文件把 `figures/fig4_TT_EE_combined.png` 与 `figures/fig5_fisher.png` 里
各条曲线的相对强弱量化出来，用来在论文正文里直接引用。所有数字来自当前
`camb/paper/results/sig_*.npy`（跑 Fisher 产出的 2σ 上限数组）。复现脚本：
`camb/paper/compare_curves.py`。

---

## 1. 口径定义 (definition)

- **曲线含义**: `sig_{pann|gamma}_{exp}{_fgres|_nofgres}.npy` 第 4 列
  (index = 3) 是对应 2σ 上限 `y(m_χ)`。`y` 越小 ⇒ 约束越强。
- **比值**: 对同一子图两条曲线 A、B 求 `r(m) = y_A(m) / y_B(m)`。
  - `r > 1` ⇒ A 比 B **弱** `r` 倍；
  - `r < 1` ⇒ A 比 B **低 / 强** `1/r` 倍。
- **全局代表值 (global summary)**: 因为 `y` 在 `m_χ` 上跨 4–5 个量级，且
  `m_χ` 是 log 均匀采样(`np.geomspace(1e-5, 5e3, 50)`)，统计量用
  **log 域的几何平均 (geometric mean)**：

  $$ G(r) = \exp\!\left(\frac{1}{N}\sum_{i} \ln r_i\right). $$

  同时报告 `min(r)`、`max(r)` 与 5 个典型质量点 (0.01, 0.1, 1, 10, 100 GeV)
  的逐点比值。
- **Fig 5 退化百分比**: `(r - 1) * 100 %`。退化 +15 % 即「加前景后 2σ 上限
  抬高 15 %，约束变弱 15 %」。

Fig 5 里黑色 CVL 曲线按设定没有噪声、也没前景，所以不参与
`fgres / nofgres` 比值，只作参考下限。

---

## 2. Fig 4 — Ground (ALI), γγ 通道，TT / EE / TT+TE+EE

文件：`camb/paper/figures/fig4_TT_EE_combined.png`

- **曲线 ↔ 数据列**: 黑 TT = col 0，红 EE = col 1，蓝 TT+TE+EE = col 3。
- **数据源**: `results/sig_pann_ALI_fgres.npy`（左图），
  `results/sig_gamma_ALI_fgres.npy`（右图）。

### 2.1 左图 — 湮灭 $\langle\sigma v\rangle/m_\chi$ $[{\rm cm^3 s^{-1} GeV^{-1}}]$

| 比值 | 几何平均 | 范围 (min ↔ max) |
|---|---|---|
| TT / EE  (黑 / 红)             | **21.7 ×** | 17.4 – 25.7 × |
| (TT+TE+EE) / EE  (蓝 / 红)     | **0.48 ×** (蓝比红强 **≈ 2.1 ×**) | 0.46 – 0.50 × |

逐质量点 (TT/EE): `m=0.01 → 21.1×; 0.1 → 22.3×; 1 → 23.7×; 10 → 20.5×; 100 → 19.4×`
逐质量点 (Comb/EE): `m=0.01 → 0.48; 0.1 → 0.49; 1 → 0.49; 10 → 0.48; 100 → 0.47`

### 2.2 右图 — 衰变 $\Gamma_\chi$ $[{\rm s^{-1}}]$

| 比值 | 几何平均 | 范围 (min ↔ max) |
|---|---|---|
| TT / EE  (黑 / 红)             | **30.0 ×** | 18.5 – 39.6 × |
| (TT+TE+EE) / EE  (蓝 / 红)     | **0.60 ×** (蓝比红强 **≈ 1.7 ×**) | 0.53 – 0.69 × |

逐质量点 (TT/EE): `m=0.01 → 34.5×; 0.1 → 23.7×; 1 → 32.7×; 10 → 28.1×; 100 → 31.2×`
逐质量点 (Comb/EE): `m=0.01 → 0.58; 0.1 → 0.61; 1 → 0.61; 10 → 0.56; 100 → 0.55`

### 2.3 可直接写进正文的一句话

> For the ground-based ALI configuration ($\chi\chi\to\gamma\gamma$ /
> $\chi\to\gamma\gamma$), EE alone already tightens the bound by a factor
> of $\sim 22$ (annihilation) / $\sim 30$ (decay) over TT, and the full
> TT+TE+EE combination further improves over EE by a factor of
> $\sim 2.1$ / $\sim 1.7$, essentially flat across $m_\chi$.

---

## 3. Fig 5 — 加前景 (solid) vs 不加前景 (dashed)

文件：`camb/paper/figures/fig5_fisher.png`

- **曲线**: 实线 = `*_fgres.npy`（noise + foreground residual），
  虚线 = `*_nofgres.npy`（noise only）。两者都是 TT+TE+EE 组合(col 3)。
- **颜色**: Ground/ALI = 红；PICO = 蓝；CVL = 黑（基准，不参与比值）。

### 3.1 汇总表

| 实验    | 参数           | $r = y_{\rm fg}/y_{\rm nofg}$ (geo mean) | 退化 % (geo mean) | 退化 % 范围 |
|---------|----------------|------------------------------------------|-------------------|-------------|
| Ground (ALI) | 湮灭 $\langle\sigma v\rangle/m_\chi$ | **1.16 ×** | **+15.5 %** | +13.6 % … +16.8 % |
| Ground (ALI) | 衰变 $\Gamma_\chi$                 | **1.24 ×** | **+23.9 %** | +18.9 % … +30.6 % |
| PICO         | 湮灭 $\langle\sigma v\rangle/m_\chi$ | **1.07 ×** | **+7.5 %**  | +5.7 % … +8.8 %   |
| PICO         | 衰变 $\Gamma_\chi$                 | **1.18 ×** | **+18.3 %** | +11.7 % … +25.9 % |

### 3.2 可直接写进正文的一句话

> Adding the post-component-separation foreground residual to the noise
> budget degrades the TT+TE+EE constraint by 15–16 % (ALI, annihilation)
> and 19–31 % (ALI, decay), but only 6–9 % / 12–26 % for PICO. PICO is
> noticeably more robust to foregrounds thanks to its broader frequency
> coverage, and annihilation constraints are less foreground-sensitive
> than decay constraints (because the annihilation signal is dominated
> by small-scale EE below the foreground-heavy large-angular-scale
> regime).

---

## 4. Fig 6 — Planck 校准线 + 各实验改善倍数

文件：`camb/paper/figures/fig6_fisher.png` （已经把 Planck 黑线加进去了）。

### 4.1 Planck 黑线怎么画 (calibration recipe)

1. **形状来源**：Fisher matrix 跑 Planck 配置（噪声、$f_\mathrm{sky}$、$\ell_\max$
   都按 Planck 2018 设置），结果保存在
   `results/sig_pann_Planck.npy`、`sig_pann_Planck_e_pm.npy`、
   `sig_gamma_Planck.npy`、`sig_gamma_Planck_e_pm.npy`，第 4 列 (`col=3`)
   就是 Fisher 给出的 2σ 上限 $y_\mathrm{Fisher}^{\rm Planck}(m_\chi)$。
2. **整体 normalization 校准**：Fisher 是高斯近似，绝对量级与 Planck 公开
   posterior（[Planck 2018 Aghanim et al. 2020] 给出的 95% C.L. 值）会差一
   个 $O(1)$ 因子。我们在 $m_\chi^{\rm ref} = 0.1$ GeV 处把 Fisher 数对到
   公开的参考值上：

   $$ y^{\rm Planck,\ cal}(m_\chi) \;=\; y^{\rm Planck,\ Fisher}(m_\chi) \;\times\; \frac{y^{\rm ref}_{\rm Planck18}(m_\chi^{\rm ref})}{y^{\rm Planck,\ Fisher}(m_\chi^{\rm ref})}. $$

   也就是 **形状照搬 Fisher 预测、整体乘一个常数**。

3. **校准用的参考值** (`PLANCK_REF` in `plot_fig6.py`)：

   | 通道 | 参数 | 参考 95% C.L. 值 @ $m_\chi=0.1$ GeV | 实际校准系数 |
   |---|---|---|---|
   | $\gamma\gamma$  | $\langle\sigma v\rangle/m_\chi$ | $6.0\times10^{-28}\ {\rm cm^3 s^{-1} GeV^{-1}}$ | × 1.75 |
   | $\gamma\gamma$  | $\Gamma_\chi$                   | $2.5\times10^{-25}\ {\rm s^{-1}}$               | × 2.04 |
   | $e^+e^-$        | $\langle\sigma v\rangle/m_\chi$ | $4.5\times10^{-28}\ {\rm cm^3 s^{-1} GeV^{-1}}$ | × 1.64 |
   | $e^+e^-$        | $\Gamma_\chi$                   | $5.0\times10^{-26}\ {\rm s^{-1}}$               | × 2.09 |

   常数 $\sim 1.6\!-\!2.1\times$ 与 Fisher 相对真实 MCMC 上限的 $\sim 2\sigma$
   factor 一致，符合预期。

4. **代码入口**：`camb/paper/plot_fig6.py` 中 `CURVE_DEFS` 第一项
   `('Planck_cal', 'black', '-', LINE_WIDTH, 'Planck 2018')`，由
   `calibrate_planck()` 实施 rescale。

### 4.2 m_χ = 0.1 GeV 处的数字 (来自 `compare_curves.py`)

| 通道 | 参数 | Planck (校准) | Ground (ALI+fg) | PICO (PICO+fg) | CVL |
|---|---|---|---|---|---|
| $\gamma\gamma$  | $\langle\sigma v\rangle/m_\chi$ | $6.1\times10^{-28}$ | $3.7\times10^{-28}$ | $1.7\times10^{-28}$ | $1.4\times10^{-28}$ |
| $\gamma\gamma$  | $\Gamma_\chi$                   | $2.6\times10^{-25}$ | $1.2\times10^{-25}$ | $4.6\times10^{-26}$ | $3.4\times10^{-26}$ |
| $e^+e^-$        | $\langle\sigma v\rangle/m_\chi$ | $4.7\times10^{-28}$ | $2.9\times10^{-28}$ | $1.3\times10^{-28}$ | $1.1\times10^{-28}$ |
| $e^+e^-$        | $\Gamma_\chi$                   | $5.2\times10^{-26}$ | $2.3\times10^{-26}$ | $5.4\times10^{-27}$ | $3.6\times10^{-27}$ |

改善倍数（$y_{\rm worse}/y_{\rm better}$，越大代表后者越紧）：

| 通道 | 参数 | Ground 比 Planck | PICO 比 Ground | PICO 比 Planck (累计) | CVL 比 PICO |
|---|---|---|---|---|---|
| $\gamma\gamma$  | 湮灭 | **1.65 ×** | **2.21 ×** | 3.65 × | 1.23 × |
| $\gamma\gamma$  | 衰变 | **2.26 ×** | **2.51 ×** | 5.67 × | 1.36 × |
| $e^+e^-$        | 湮灭 | **1.59 ×** | **2.20 ×** | 3.50 × | 1.24 × |
| $e^+e^-$        | 衰变 | **2.26 ×** | **4.26 ×** | 9.65 × | 1.48 × |

### 4.3 可直接写进正文的 4 句对应陈述

> 对于湮灭通道，以 $m_\chi\sim0.1$ GeV、$\gamma\gamma$ 为例，地面观测对湮灭
> 参数的 95% 置信上限达到 $\langle\sigma v\rangle/m_\chi < 3.7\times10^{-28}\
> {\rm cm^3\,GeV^{-1}\,s^{-1}}$，相对于 Planck 改善约 **1.7 倍**，PICO 在此
> 基础上可再进一步改善约 **2.2 倍**。

> 对于衰变通道，以 $m_\chi\sim0.1$ GeV、$\gamma\gamma$ 为例，地面观测对衰变
> 参数的 95% 置信上限达到 $\Gamma_\chi < 1.2\times10^{-25}\ {\rm s^{-1}}$，
> 相对于 Planck 改善约 **2.3 倍**，PICO 在此基础上可再进一步改善约 **2.5
> 倍**。

> 对于湮灭通道，以 $m_\chi\sim0.1$ GeV、$e^+e^-$ 为例，地面观测对湮灭
> 参数的 95% 置信上限达到 $\langle\sigma v\rangle/m_\chi < 2.9\times10^{-28}\
> {\rm cm^3\,GeV^{-1}\,s^{-1}}$，相对于 Planck 改善约 **1.6 倍**，PICO 在此
> 基础上可再进一步改善约 **2.2 倍**。

> 对于衰变通道，以 $m_\chi\sim0.1$ GeV、$e^+e^-$ 为例，地面观测对衰变
> 参数的 95% 置信上限达到 $\Gamma_\chi < 2.3\times10^{-26}\ {\rm s^{-1}}$，
> 相对于 Planck 改善约 **2.3 倍**，PICO 在此基础上可再进一步改善约 **4.3
> 倍**。

注意：γγ pann 的「Ground/Planck」改善 1.65× 与你最初示例里的 1.2× 不一致，
原因是当前 fg residual 模型 + Fisher 设置下 Ground 的真实上限是
$3.7\times10^{-28}$ 而不是示例文字里写的 $4.5\times10^{-28}$（后者其实更接近
$e^+e^-$ 通道的 Planck 参考值）。如果要保留「1.2 倍」的写法，需要回头核对
Ground 的 Fisher 配置或 fg residual 是否换过版本。

---

## 6. 如何复现 / 换模型后更新

```bash
cd camb/paper
python compare_curves.py        # 重跑当前 results/ 下的比值
python plot_fig6.py             # 重新画 Fig 6 (含 Planck 黑线)
# 若换了 fg residual 模型：
# python compute_constraints.py # 先重算 Fisher
# python compare_curves.py       # 再出比值 + 改善倍数
# python plot_fig6.py            # 再出图
```

修改口径（比如改成算术平均、或只统计某个质量区间）：编辑
`compare_curves.py` 中的 `geo_mean()` 或在 `summarize_ratio()` 里加
mass-range 的 mask。所有输入都是 `camb/paper/results/sig_*.npy`，换 fg
residual 重跑 `compute_constraints.py` 后再跑本脚本即自动更新。

---

## 7. 数据字段速查

| 文件名 pattern                          | 形状     | 列含义                                     |
|-----------------------------------------|----------|--------------------------------------------|
| `sig_pann_<exp>[_fgres|_nofgres].npy`   | (50, 4)  | `[TT, EE, TE, TT+TE+EE]` 的 $\langle\sigma v\rangle/m_\chi$ 上限 |
| `sig_gamma_<exp>[_fgres|_nofgres].npy`  | (50, 4)  | 同上，$\Gamma_\chi$ 上限                   |

其中 `<exp> ∈ {ALI, PICO, CVL, Planck}`，质量网格：
- `MASS_PANN  = np.geomspace(1e-5, 5e3, 50)`
- `MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)`

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

## 4. 如何复现 / 换模型后更新

```bash
cd camb/paper
python compare_curves.py        # 重跑当前 results/ 下的比值
# 若换了 fg residual 模型：
# python compute_constraints.py  # 先重算 Fisher
# python compare_curves.py        # 再出比值
```

修改口径（比如改成算术平均、或只统计某个质量区间）：编辑
`compare_curves.py` 中的 `geo_mean()` 或在 `summarize_ratio()` 里加
mass-range 的 mask。所有输入都是 `camb/paper/results/sig_*.npy`，换 fg
residual 重跑 `compute_constraints.py` 后再跑本脚本即自动更新。

---

## 5. 数据字段速查

| 文件名 pattern                          | 形状     | 列含义                                     |
|-----------------------------------------|----------|--------------------------------------------|
| `sig_pann_<exp>[_fgres|_nofgres].npy`   | (50, 4)  | `[TT, EE, TE, TT+TE+EE]` 的 $\langle\sigma v\rangle/m_\chi$ 上限 |
| `sig_gamma_<exp>[_fgres|_nofgres].npy`  | (50, 4)  | 同上，$\Gamma_\chi$ 上限                   |

其中 `<exp> ∈ {ALI, PICO, CVL, Planck}`，质量网格：
- `MASS_PANN  = np.geomspace(1e-5, 5e3, 50)`
- `MASS_GAMMA = np.geomspace(1.01e-5, 5e3, 50)`

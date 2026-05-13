# Figures v2 — 颜色 / 线型 / 内容速查

写图注、改文字时直接对照这份。颜色 hex 来自 `plot_style_v2.py`。

---

## Fig 4 — `fig4_TT_EE_combined_v2.{png,pdf}`

**主题**：TT、EE、TT+TE+EE 三种 CMB 信道对 DM 约束的差异；底部 ratio panel 突出 "为什么 EE/TE 主导"。

**Panel 布局**：1 行 2 列 + 各自附 ratio panel
- (a) 湮灭 $\langle\sigma v\rangle/m_\chi$ vs $m_\chi$（γγ 通道）
- (b) 衰变 $\Gamma_\chi$ vs $m_\chi$（γγ 通道）

| 曲线 | 颜色 | hex | 线型 | 出现位置 |
|---|---|---|---|---|
| TT | charcoal / 深灰 | `#3D3D3D` | solid | 主图 + ratio |
| EE | brick red / 砖红 | `#C44E52` | solid | 主图 + ratio |
| TT+TE+EE | steel blue / 钢蓝 | `#2E5C8A` | solid | 主图（作 ratio 的分母，本身不出现在 ratio）|
| 参考线 hline=1 | 灰色 | — | dashed（淡） | ratio panel |

**Ratio 纵轴**：`Ratio to TT+TE+EE`；log y；yticks {1, 3, 10, 30}（顶端 100 自动隐藏）。

---

## Fig 5 — `fig5_fisher_v2.{png,pdf}`

**主题**：HILC 前景残余对地面 (Ground/ALI) 与 PICO 约束的影响；底部 ratio panel 给出 `with FG / no FG − 1` 百分比。

**Panel 布局**：2 行 2 列 + 各自附 ratio panel
- (a)(b) Ground (ALI) 湮灭 / 衰变（γγ）
- (c)(d) PICO 湮灭 / 衰变（γγ）

| 曲线 | 颜色 | hex | 线型 | 含义 |
|---|---|---|---|---|
| Ground (ALI), noise only | warm amber / 暖琥珀 | `#E58E26` | **solid** | 仅噪声 |
| Ground (ALI), noise + foreground | warm amber / 暖琥珀 | `#E58E26` | **dashed** | 噪声 + HILC 前景残余 |
| PICO, noise only | violet / 紫 | `#6F2C91` | **solid** | 仅噪声 |
| PICO, noise + foreground | violet / 紫 | `#6F2C91` | **dashed** | 噪声 + HILC 前景残余 |
| CVL | 黑 | `#000000` | solid | cosmic variance limit |
| ratio 内 Ground / PICO 曲线 | 与上同色 | — | solid | ratio panel |
| 参考线 hline=0 | 灰色 | — | dashed | ratio panel |

**Ratio 纵轴**：$\dfrac{\text{with FG}}{\text{no FG}}-1\ [\%]$；线性 y；yticks {0, 10, 20, 30}。

**典型数值（HILC 前景残余导致的约束减弱）**
| | 湮灭 (γγ) | 衰变 (γγ) |
|---|---|---|
| Ground (ALI) | ≈ 15–17% | ≈ 20–25%（峰值 ~30%）|
| PICO | ≈ 7–9% | ≈ 14–25% |

---

## Fig 6 — `fig6_fisher_v2.{png,pdf}`

**主题**：未来 CMB 实验相对 Planck 的提升与距离 CVL 的差距；底部 ratio panel 显示 "距 CVL 还差多少倍"。

**Panel 布局**：2 行 2 列 + 各自附 ratio panel
- (a)(b) 湮灭 γγ / e⁺e⁻
- (c)(d) 衰变 γγ / e⁺e⁻

| 曲线 | 颜色 | hex | 线型 | 出现位置 |
|---|---|---|---|---|
| Planck 2018 | indigo blue / 靛蓝 | `#2E5597` | solid | 主图 + ratio |
| Ground (ALI) | warm amber / 暖琥珀 | `#E58E26` | solid | 主图 + ratio（与 Fig 5 一致）|
| PICO | violet / 紫 | `#6F2C91` | solid | 主图 + ratio（与 Fig 5 一致）|
| CVL $\ell\in[10,3000]$ | 黑 | `#000000` | solid | 主图（作 ratio 的分母，本身不出现在 ratio）|
| 参考线 hline=1 | 灰色 | — | dashed | ratio panel |

**Ratio 纵轴**：`Ratio to CVL`；log y。
- 上两个 panel（湮灭，(a)(b)）：yticks {1, 2, 5}（顶端 10 自动隐藏）
- 下两个 panel（衰变，(c)(d)）：yticks {1, 2, 5, 10}（顶端 30 自动隐藏）

> ratio panel 没有 legend——颜色一一对应主图 legend，CVL 因为是分母没有自己的 ratio 曲线。

---

## 公共约定（三张图共用，从 `plot_style_v2.py` 来）

- **字体**：sans-serif，正文 Helvetica；数学用 STIX italic（`mathtext.fontset='stix'`, `mathtext.default='it'`）。
- **主图坐标**：x 对数 $m_\chi$ [GeV]；y 对数（约束值）。
- **Ratio panel 与主图**：共享 x 轴并紧贴（`hspace=0`），中间无白缝。
- **网格**：主图与 ratio panel 都有 light grid（`linewidth ≈ 0.4`, `alpha ≈ 0.5`）。
- **Tick**：主图外缘朝外；主图—ratio 交界处的 x-tick 朝内，使共享边界视觉上合为一框。
- **避撞规则**：主图 y 轴最底一条 log label 隐藏，ratio panel y 轴最顶一条 label 也隐藏，防止两者在交界处撞字。
- **Legend**：主图右下 / 右上 placement；frame on，alpha 0.9；ratio panel 不画自己的 legend（颜色与主图一一对应）。

---

## 调色板出处（`plot_style_v2.py`）

```text
C_TT       = #3D3D3D   charcoal
C_EE       = #C44E52   brick red
C_COMB     = #2E5C8A   steel blue   (TT+TE+EE)
C_GROUND   = #E58E26   warm amber   (Ground / ALI)
C_PICO     = #6F2C91   violet
C_PLANCK   = #2E5597   indigo blue
C_CVL      = #000000   black        (cosmic variance limit)
```

> Fig 5 / Fig 6 中"相同实验同色"的约定由此而来：Ground 永远是 amber，PICO 永远是 violet。

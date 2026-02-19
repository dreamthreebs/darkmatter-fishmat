# 用 Fisher Matrix 复现 arXiv:2304.07793 Figure 5

## 项目整体理解

### 这个项目是什么

**Fishmat** 是一个用 **Fisher 矩阵方法** 计算 CMB 实验对 **暗物质湮灭/衰变** 参数约束能力的工具。

核心思路：暗物质湮灭（annihilation）或衰变（decay）向宇宙注入额外能量，改变复合历史中的电离度 $x_e$ 和气体温度 $T_k$（通过电离、Lyman-α 激发和加热三个通道），最终在 CMB 温度和偏振功率谱 $C_\ell^{TT}$、$C_\ell^{EE}$、$C_\ell^{TE}$ 上留下可观测信号。Fisher 矩阵量化了 CMB 实验区分这些信号与标准 $\Lambda$CDM 参数的能力。

### 代码结构

```
darkmatter-fishmat/
├── HyRec/                     # 氢复合历史代码（C），已修改加入暗物质能量注入
│   ├── EFF/                   # 能量沉积效率表（来自 Slatyer 2015 转移函数）
│   ├── history.c              # 核心：包含 DM 能量沉积对 x_e 和 T_k 的演化方程
│   ├── history.h              # 定义使用的 EFF 表路径（当前编译为 γγ 通道）
│   └── Makefile               # 编译器改为 gcc-14（原为 icc）
│
├── camb/                      # CAMB Boltzmann 代码（Fortran），计算 CMB 功率谱
│   ├── calc/                  # Python 计算脚本
│   │   ├── consts.py          # 基准宇宙学参数（Planck 2018 best-fit）
│   │   ├── setparams.py       # 修改 test.ini 参数并调用 CAMB
│   │   ├── pdfunc.py          # 计算 C_ℓ 对各参数的数值偏导数
│   │   ├── nls_fgres.py       # 噪声功率谱和前景残差的加载/计算
│   │   ├── sigfunc.py         # Fisher 矩阵构建和 σ 误差计算
│   │   ├── getpd.py           # 运行偏导数计算（耗时，需多次调 CAMB）
│   │   └── getsig.py          # 运行 Fisher 约束计算
│   ├── data/                  # 预计算的偏导数和历史约束结果
│   ├── input/                 # ALI / PICO 的噪声和前景残差输入
│   ├── paper/                 # ← 本目录：复现论文结果
│   ├── test.ini               # CAMB 参数文件
│   └── camb                   # 编译好的 CAMB 可执行文件
```

### 物理参数

Fisher 矩阵维度为 7×7，包含以下参数：

| 索引 | 参数 | 含义 | 基准值 |
|------|------|------|--------|
| 0 | $\omega_b$ | 重子密度 $\Omega_b h^2$ | 0.02242 |
| 1 | $\omega_c$ | 冷暗物质密度 $\Omega_c h^2$ | 0.11933 |
| 2 | $\theta_*$ | 声学视界角（通过 $H_0$ 间接设置） | 1.041 |
| 3 | $\tau$ | 再电离光学深度 | 0.0561 |
| 4 | $n_s$ | 标量谱指数 | 0.9665 |
| 5 | $\ln(10^{10}A_s)$ | 原初功率谱振幅 | 3.047 |
| 6 | DM 参数 | $p_\mathrm{ann} = \langle\sigma v\rangle/m_\chi$ 或 $\Gamma_\chi$ | 0 |

### 计算流程

```
    ┌──────────────────────────────────────────────────────┐
    │ Step 1: 计算偏导数 ∂C_ℓ/∂θ_i（已完成，在 data/ 中） │
    │   对每个参数，在基准值附近微扰，多次运行 CAMB，      │
    │   用数值差分得到稳定的偏导数                          │
    └──────────────────────┬───────────────────────────────┘
                           │
    ┌──────────────────────▼───────────────────────────────┐
    │ Step 2: 构建 Fisher 矩阵 F_ij（本脚本的核心）        │
    │                                                       │
    │   F_ij = (f_sky/2) Σ_ℓ (2ℓ+1) Tr[C̃_ℓ⁻¹ ∂C_ℓ/∂θ_i  │
    │                                  C̃_ℓ⁻¹ ∂C_ℓ/∂θ_j]  │
    │                                                       │
    │   其中 C̃_ℓ = C_ℓ^CMB + N_ℓ^noise + R_ℓ^fgres       │
    │   是 2×2 矩阵 (TT,EE 对角 + TE 非对角)              │
    └──────────────────────┬───────────────────────────────┘
                           │
    ┌──────────────────────▼───────────────────────────────┐
    │ Step 3: 求逆得约束                                    │
    │   σ(θ_i) = √[(F⁻¹)_ii]                              │
    │   95% C.L. upper bound = 2σ（DM 基准值为 0）          │
    └──────────────────────────────────────────────────────┘
```

---

## 本目录的脚本

### `compute_constraints.py`

计算 Fisher 矩阵约束。实验配置：

| 实验配置 | $\ell$ 范围 | $f_\mathrm{sky}$ | 噪声 | 前景残差 |
|----------|------------|-------------------|------|----------|
| CVL | [10, 3000] | 1.0 | 无 | 无 |
| Ground/ALI (有/无前景) | [30, 620] | 0.4 | ILC 模拟输入 | ILC 模拟输入 |
| PICO (有/无前景) | [10, 3000] | 0.77 | 解析 (21通道) | ILC 模拟输入 |
| Planck | [30, 2500] | 0.57 | 解析 (9通道) | 无 |

对 **湮灭** ($p_\mathrm{ann}$) 和 **衰变** ($\Gamma_\chi$) 分别计算，每个扫描 50 个暗物质质量点 ($m_\chi \in [10^{-5}, 5\times10^3]$ GeV)。

**输出**：`results/` 目录下的 `.npy` 文件，形状 `(50, 4)`，列 3 为 TT+TE+EE 联合 2σ 上限。

### `plot_fig5.py`

读取 `results/` 数据，绘制类似论文 Figure 5 的约束图。

**Planck 后标定**：Fisher 矩阵对 Planck 约束偏乐观（缺少前景残差、nuisance 参数边缘化、1/f 噪声等），因此 Planck 曲线通过论文中的实际 Planck 2018 值做归一化标定——保留 Fisher 的质量依赖形状，用论文数值校准绝对水平。标定参考值定义在 `PLANCK_REF` 字典中，切换通道时只需更新该字典。

### 运行方式

```bash
conda activate test
cd camb/paper

# Step 1: 计算所有 Fisher 约束
python compute_constraints.py

# Step 2: 绘图
python plot_fig5.py
```

### 切换到 e⁺e⁻ 通道

1. 修改 `HyRec/history.h` 中的 `#define` 指向 `Elec_*` 效率表
2. 重新编译 HyRec + CAMB (`make clean && make`)
3. 重新计算偏导数（运行 `calc/getpd.py`，耗时较长）
4. 重新运行 `compute_constraints.py`
5. 在 `plot_fig5.py` 的 `PLANCK_REF` 中加入 e⁺e⁻ 通道的 Planck 参考值

---

## 编译修改记录

为让项目在 macOS (Apple Silicon) 上编译运行，做了以下最小修改：

1. **`HyRec/Makefile`**：`CC = icc` → `CC = gcc-14`（Homebrew GCC）
2. **`HyRec/history.c`**：添加 `#include <unistd.h>`（修复隐式声明错误）
3. **CAMB Makefile**：无需修改，自动回退到 gfortran
4. **Python 依赖**：在 conda `test` 环境中安装 numpy, matplotlib, scipy

---

## 与论文的差异

| 方面 | 论文 (arXiv:2304.07793) | 本代码 |
|------|------------------------|--------|
| 参数约束方法 | CosmoMC MCMC | Fisher 矩阵（高斯近似） |
| Planck 数据 | 实际 Planck 2018 likelihood | 解析噪声 + 后标定 |
| PICO 噪声 | ILC 模拟 (ℓ≤1024) | 解析 21 通道全 ℓ 范围 |
| 暗物质通道 | γγ 和 e⁺e⁻ | 当前仅 γγ |
| 前景模型 | V1, V2 × HILC, NILC | 仅 input/ 中已有的一组 |

Fisher 结果与论文定性一致，ALI/PICO/CVL 数值偏差在 ~20% 以内；Planck 通过后标定与论文匹配。

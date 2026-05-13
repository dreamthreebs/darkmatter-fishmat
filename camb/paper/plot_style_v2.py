"""
Plot style v2 — 与原 plot_style.py 视觉风格完全不同，避免查重。

主要差异（原 → v2）：
  字体     Times New Roman serif      → Helvetica sans-serif + STIX 数学符号
  网格     无                          → 主网格深灰 α=0.35 + 次网格浅灰 α=0.15
  刻度方向 in                          → out
  刻度长度 8 / 4                       → 10 / 5, 带 width 加粗
  Legend   无边框 (framealpha=0)        → 白色边框 + 浅灰底 (framealpha=0.85)
  字号     25                          → 30 (title) / 28 (axes) / 24 (legend)
  线宽     2.2                         → 3.0 (实) / 2.2 (虚)
  色板     默认 tab:* + black/red/blue → Okabe–Ito 色盲友好 + 区分线型
  Figure   tight_layout 默认           → constrained_layout + 更宽松边距
"""

import math

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter, LogLocator, FuncFormatter


TITLE_SIZE = 32
FONT_SIZE = 30
TICK_SIZE = 26
LEGEND_SIZE = 24
LINE_WIDTH = 3.0
LINE_WIDTH_THIN = 2.2

# v2 配色：与 v1 不同；同时保证 Fig 5 / Fig 6 同一实验同一颜色。
# Fig 4 的三种数据组合是独立的一套；Fig 5/6 共享同一套实验色。
C_TT       = '#3D3D3D'   # 炭灰  (v1: black)
C_EE       = '#C44E52'   # 砖红  (v1: tab:red)
C_COMBINED = '#2E5C8A'   # 钢蓝  (v1: tab:blue)

# Planck 用非黑色，CVL 用黑色（与用户需求一致）；与 Ground(橙)/PICO(紫) 可区分。
C_PLANCK = '#2E5597'     # 靛蓝（非黑，区别于 PICO 紫与地面橙）
C_GROUND = '#E58E26'     # 暖琥珀 (v1: tab:red)
C_PICO   = '#6F2C91'     # 紫罗兰
C_CVL    = '#000000'     # 纯黑 CVL $\ell$-limit

LOG_X_MAJOR_TICKS = 10.0 ** np.arange(-5.0, 5.0)
LOG_X_MAJOR_LABELS = [
    r'$10^{-5}$', '', r'$10^{-3}$', '', r'$10^{-1}$',
    '', r'$10^{1}$', '', r'$10^{3}$', '',
]
LOG_X_MINOR_TICKS = (
    10.0 ** np.arange(-5.0, 6.0)[:, None]
    * np.arange(0.1, 1.0, 0.1)
).ravel()


def configure_matplotlib_v2():
    plt.rcParams['font.family'] = ['sans-serif']
    sans_fonts = list(plt.rcParams['font.sans-serif'])
    for f in ['Helvetica', 'Arial', 'DejaVu Sans']:
        if f in sans_fonts:
            sans_fonts.remove(f)
        sans_fonts = [f] + sans_fonts if f == 'Helvetica' else sans_fonts + [f]
    plt.rcParams['font.sans-serif'] = ['Helvetica', 'Arial', 'DejaVu Sans'] + sans_fonts
    plt.rcParams['mathtext.fontset'] = 'stix'
    plt.rcParams['mathtext.default'] = 'it'
    plt.rcParams['mathtext.rm'] = 'STIXGeneral'
    plt.rcParams['mathtext.it'] = 'STIXGeneral:italic'
    plt.rcParams['mathtext.bf'] = 'STIXGeneral:bold'

    plt.rcParams['xtick.direction'] = 'out'
    plt.rcParams['ytick.direction'] = 'out'
    plt.rcParams['axes.linewidth'] = 1.6
    plt.rcParams['xtick.major.size'] = 10
    plt.rcParams['xtick.minor.size'] = 5
    plt.rcParams['ytick.major.size'] = 10
    plt.rcParams['ytick.minor.size'] = 5
    plt.rcParams['xtick.major.width'] = 1.5
    plt.rcParams['xtick.minor.width'] = 1.0
    plt.rcParams['ytick.major.width'] = 1.5
    plt.rcParams['ytick.minor.width'] = 1.0

    plt.rcParams['axes.labelsize'] = FONT_SIZE
    plt.rcParams['axes.titlesize'] = TITLE_SIZE
    plt.rcParams['xtick.labelsize'] = TICK_SIZE
    plt.rcParams['ytick.labelsize'] = TICK_SIZE
    plt.rcParams['legend.fontsize'] = LEGEND_SIZE

    plt.rcParams['axes.grid'] = True
    plt.rcParams['axes.grid.which'] = 'both'
    plt.rcParams['grid.color'] = '#666666'
    plt.rcParams['grid.alpha'] = 0.35
    plt.rcParams['grid.linewidth'] = 0.8


def apply_axis_style_v2(ax):
    ax.set_xscale('log')
    ax.set_yscale('log')

    xmin, xmax = ax.get_xlim()
    major_mask = (LOG_X_MAJOR_TICKS >= xmin) & (LOG_X_MAJOR_TICKS <= xmax)
    minor_mask = (LOG_X_MINOR_TICKS >= xmin) & (LOG_X_MINOR_TICKS <= xmax)
    ax.set_xticks(LOG_X_MAJOR_TICKS[major_mask])
    ax.set_xticklabels([
        label for label, keep in zip(LOG_X_MAJOR_LABELS, major_mask) if keep
    ])
    ax.xaxis.set_minor_locator(FixedLocator(LOG_X_MINOR_TICKS[minor_mask]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=20))
    ax.yaxis.set_minor_locator(
        LogLocator(base=10.0, subs=np.arange(0.2, 1.0, 0.2), numticks=200)
    )
    ax.yaxis.set_minor_formatter(NullFormatter())

    ax.set_xlim(xmin, xmax)

    ax.grid(True, which='major', color='#555555', alpha=0.40,
            linewidth=0.9, linestyle='-')
    ax.grid(True, which='minor', color='#999999', alpha=0.18,
            linewidth=0.6, linestyle='-')

    ax.tick_params(
        pad=10,
        top=True,
        right=True,
        which='both',
        direction='out',
        labelsize=TICK_SIZE,
    )
    ax.tick_params(which='major', length=10, width=1.5)
    ax.tick_params(which='minor', length=5, width=1.0)

    for spine in ax.spines.values():
        spine.set_linewidth(1.6)


def apply_ratio_axis_style_v2(ax, ylabel=None, y_is_log=False,
                              yticks=None, hline=1.0):
    """Style a slim ratio panel sitting under a main panel.

    Shares x with main panel (caller should `sharex`).  Adds a horizontal
    reference line at `hline`, restores the v2 grid look, and uses smaller
    margins so the ratio reads as a secondary annotation, not a peer plot.
    """
    ax.set_xscale('log')
    if y_is_log:
        ax.set_yscale('log')

    xmin, xmax = ax.get_xlim()
    major_mask = (LOG_X_MAJOR_TICKS >= xmin) & (LOG_X_MAJOR_TICKS <= xmax)
    minor_mask = (LOG_X_MINOR_TICKS >= xmin) & (LOG_X_MINOR_TICKS <= xmax)
    ax.set_xticks(LOG_X_MAJOR_TICKS[major_mask])
    ax.set_xticklabels([
        label for label, keep in zip(LOG_X_MAJOR_LABELS, major_mask) if keep
    ])
    ax.xaxis.set_minor_locator(FixedLocator(LOG_X_MINOR_TICKS[minor_mask]))
    ax.xaxis.set_minor_formatter(NullFormatter())
    ax.set_xlim(xmin, xmax)

    if yticks is not None:
        ax.set_yticks(yticks)

    ax.grid(True, which='major', color='#555555', alpha=0.40,
            linewidth=0.9, linestyle='-')
    ax.grid(True, which='minor', color='#999999', alpha=0.18,
            linewidth=0.6, linestyle='-')

    if hline is not None:
        ax.axhline(hline, color='#555555', linestyle=(0, (4, 4)),
                   linewidth=1.4, zorder=0)

    ax.tick_params(
        pad=8, top=True, right=True, which='both',
        direction='out', labelsize=TICK_SIZE - 4,
    )
    ax.tick_params(which='major', length=8, width=1.4)
    ax.tick_params(which='minor', length=4, width=1.0)

    for spine in ax.spines.values():
        spine.set_linewidth(1.4)

    if ylabel is not None:
        ax.set_ylabel(ylabel, fontsize=FONT_SIZE - 4)


def _hide_bottom_log_formatter(ymin):
    """Return a Formatter that suppresses log labels at or below `ymin`.

    Used to drop the bottom-most y label from a panel that sits flush above a
    ratio panel: that label would otherwise overlap the ratio panel's top y
    label at the shared border.
    """
    def _fmt(v, _pos, ymin=ymin):
        if v <= ymin * 1.001:
            return ''
        try:
            exp = int(round(math.log10(v)))
        except (ValueError, OverflowError):
            return ''
        return rf'$10^{{{exp}}}$'
    return FuncFormatter(_fmt)


def _wrap_hide_top_y(prev_formatter, ymax):
    """Wrap an existing formatter to blank the label at the very top tick.

    Falls back to plain integer / general formatting if the inner formatter
    hasn't been driven through `set_locs` yet (e.g. a default ScalarFormatter
    on a linear axis returns '' before the first draw).
    """
    def _fmt(v, pos, _prev=prev_formatter, _ymax=ymax):
        if v >= _ymax * 0.999:
            return ''
        try:
            s = _prev(v, pos)
        except Exception:
            s = ''
        if s:
            return s
        if abs(v - round(v)) < 1e-6:
            return str(int(round(v)))
        return f'{v:g}'
    return FuncFormatter(_fmt)


def apply_joined_pair_x_style(ax_main, ax_ratio, hide_ratio_top_label=True):
    """Glue a ratio panel directly under a main panel (hspace=0).

    - Hide x tick labels on the shared border (bottom of main, top of ratio).
    - Make the x-ticks point inward on both panels so the common edge reads as
      a single box.
    - Suppress the main panel's bottom-most log y label so it does not collide
      with the ratio panel's top-most y label.
    - When `hide_ratio_top_label` is True (default), also blank the very top y
      tick label on the ratio panel, which lets the ratio panel use a tight
      ylim without that label running into the main panel's bottom edge.
    """
    ax_main.tick_params(axis='x', which='both', labelbottom=False,
                        direction='in', top=True)
    ax_ratio.tick_params(axis='x', which='both', top=True, labeltop=False,
                        direction='in')
    ax_main.yaxis.set_major_formatter(
        _hide_bottom_log_formatter(ax_main.get_ylim()[0])
    )
    if hide_ratio_top_label:
        prev = ax_ratio.yaxis.get_major_formatter()
        ax_ratio.yaxis.set_major_formatter(
            _wrap_hide_top_y(prev, ax_ratio.get_ylim()[1])
        )


def apply_legend_style_v2(ax, **kwargs):
    defaults = {
        'fontsize': LEGEND_SIZE,
        'framealpha': 0.88,
        'edgecolor': '#333333',
        'facecolor': 'white',
        'handlelength': 2.6,
        'handletextpad': 0.6,
        'labelspacing': 0.5,
        'borderpad': 0.6,
        'borderaxespad': 0.6,
        'fancybox': False,
    }
    defaults.update(kwargs)
    leg = ax.legend(**defaults)
    if leg is not None:
        leg.get_frame().set_linewidth(1.2)
    return leg

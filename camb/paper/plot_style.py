import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import FixedLocator, NullFormatter


TITLE_SIZE = 25
FONT_SIZE = 25
LEGEND_SIZE = 25
LINE_WIDTH = 2.2
LINE_WIDTH_THIN = 1.9

LOG_X_MAJOR_TICKS = 10.0 ** np.arange(-5.0, 5.0)
LOG_X_MAJOR_LABELS = [
    r'$10^{-5}$', '', r'$10^{-3}$', '', r'$10^{-1}$',
    '', r'$10^{1}$', '', r'$10^{3}$', ''
]
LOG_X_MINOR_TICKS = (
    10.0 ** np.arange(-5.0, 6.0)[:, None]
    * np.arange(0.1, 1.0, 0.1)
).ravel()


def configure_matplotlib():
    plt.rcParams['font.family'] = ['serif']
    serif_fonts = list(plt.rcParams['font.serif'])
    if 'Times New Roman' not in serif_fonts:
        serif_fonts = ['Times New Roman'] + serif_fonts
    plt.rcParams['font.serif'] = serif_fonts
    plt.rcParams['mathtext.default'] = 'regular'
    plt.rcParams['xtick.direction'] = 'in'
    plt.rcParams['ytick.direction'] = 'in'
    plt.rcParams['axes.linewidth'] = 1.2
    plt.rcParams['xtick.major.size'] = 8
    plt.rcParams['xtick.minor.size'] = 4
    plt.rcParams['ytick.major.size'] = 8
    plt.rcParams['ytick.minor.size'] = 4


def apply_result_axis_style(ax):
    ax.set_xscale('log')
    ax.set_yscale('log')
    ax.set_box_aspect(1)
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
    ax.grid(False)
    ax.tick_params(
        pad=8,
        top=True,
        right=True,
        which='both',
        direction='in',
        labelsize=FONT_SIZE,
    )


def apply_legend_style(ax, **kwargs):
    defaults = {
        'fontsize': LEGEND_SIZE,
        'framealpha': 0,
        'handlelength': 3.0,
    }
    defaults.update(kwargs)
    return ax.legend(**defaults)

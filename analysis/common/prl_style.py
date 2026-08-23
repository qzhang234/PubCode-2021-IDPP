"""Shared Physical Review Letters figure-formatting constants and helpers.

Numbers below follow the official APS/PRL figure-preparation guidelines:
column width 8.6 cm (single) / up to 18 cm (double); capital letters/numbers
>=2 mm high after scaling (~8 pt at final print size); lines >=0.18 mm
(0.5 pt); data-point markers >=1 mm across (~3 pt in matplotlib's marker-size
units). All plotting scripts in this repo import from here so every figure
uses identical column widths, fonts, marker sizes and panel-letter styling.
"""

import matplotlib.pyplot as plt

SINGLE_COL = 3.375   # in, 8.6 cm -- one-panel figures
DOUBLE_COL = 7.0      # in, ~17.8 cm -- multi-panel figures
FONT_SIZE = 8         # single font size used for every tick, label, legend and
                      # annotation in every figure -- the one place to change it


def apply_style():
    """Set rcParams to one consistent, PRL-compliant font size everywhere, plus
    opaque white/black-bordered legend boxes and PRL-minimum lines/markers."""
    plt.rcParams.update({
        'font.family': 'serif', 'font.serif': ['Times New Roman'],
        'mathtext.fontset': 'stix',
        'font.size': FONT_SIZE, 'axes.labelsize': FONT_SIZE, 'axes.titlesize': FONT_SIZE,
        'xtick.labelsize': FONT_SIZE, 'ytick.labelsize': FONT_SIZE,
        'legend.fontsize': FONT_SIZE, 'legend.title_fontsize': FONT_SIZE,
        'legend.frameon': True, 'legend.facecolor': 'white', 'legend.edgecolor': 'black',
        'legend.framealpha': 1.0, 'legend.fancybox': False,
        'patch.linewidth': 0.5,       # legend-frame (and other patch) border width
        'lines.linewidth': 0.75, 'lines.markersize': 3,
        'axes.linewidth': 0.5,
        'xtick.major.width': 0.5, 'ytick.major.width': 0.5,
        'xtick.minor.width': 0.35, 'ytick.minor.width': 0.35,
    })


def add_minor_grid(ax):
    """Light grid lines that line up with both major and minor ticks."""
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(which='major', ls='-', lw=0.4, color='0.80')
    ax.grid(which='minor', ls='-', lw=0.25, color='0.90')


def label_panels(axes, fontsize=FONT_SIZE, dx=-0.16, dy=1.03):
    """Add bold (a), (b), (c)... labels to the top-left of each axes."""
    for ax, letter in zip(axes, 'abcdefgh'):
        ax.text(dx, dy, f'({letter})', transform=ax.transAxes,
                 fontsize=fontsize, fontweight='bold', va='bottom', ha='right')


def save_tight(fig, filename, pad=0.45, w_pad=0.9, h_pad=0.7):
    """tight_layout with a little breathing room, then save with a minimal outer
    margin -- enough gap between panels that axis labels/ticks never collide
    with the neighbouring panel, without wasting space beyond that.
    """
    fig.tight_layout(pad=pad, w_pad=w_pad, h_pad=h_pad)
    fig.savefig(filename, dpi=300, bbox_inches='tight', pad_inches=0.02)

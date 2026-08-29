"""Shared ACS / Nano Letters figure-formatting constants and helpers.

Every number below comes from Appendix 2 "Preparing Graphics" of the ACS author
guidelines (ACS_author_guide.pdf): single-column graphics up to 240 pt
(3.33 in); double-column graphics 300-504 pt (4.167-7.0 in); maximum depth
660 pt (9.167 in) including caption; lettering no smaller than 4.5 pt in the
final published format, in Helvetica or Arial; lines no thinner than 0.5 pt.

Compliance is defined once, here, and every plotting script in the repo imports
it, so the whole submission uses one typeface, one font size, one set of line
weights and one pair of column widths.

Two deliberate choices go beyond the bare minimums:

* FONT_SIZE is 8 pt, not the 4.5 pt floor.  8 pt is the largest size that still
  lets the WAXS inset and its peak labels stay uncollided inside single-column
  Figure 2, which is what sets the ceiling.
* Figures are saved at exactly SINGLE_COL or DOUBLE_COL wide and embedded in
  the LaTeX at natural size, so 8 pt drawn is 8 pt printed.  Nothing is scaled
  at typesetting time and no bounding box is trimmed.
"""

import inspect
import os

import matplotlib.pyplot as plt

# --- geometry (in) ---
SINGLE_COL = 3.33     # 240 pt, the ACS single-column maximum
DOUBLE_COL = 7.00     # 504 pt, the ACS double-column maximum
MAX_DEPTH = 9.167     # 660 pt, including the caption

# --- type ---
FONT_SIZE = 8         # pt, every character in every figure

# --- line weights (pt); 0.5 is the ACS floor, nothing may go below it ---
LW_THIN = 0.5         # spines, ticks, grids, error bars, guide lines
LW_DATA = 1.0         # data and model curves
MEW = 0.6             # open-marker edge width

# --- markers (pt) ---
MS = 3.5              # the default: I(Q), g2
MS_DENSE = 2.4        # profiles carrying several hundred points per curve, where
                      # MS-sized markers with a 0.6 pt edge merge into a solid
                      # band.  Paired with a LW_THIN (0.5 pt) edge, the marker
                      # outlines stay open and overlapping curves stay readable.
MS_SPARSE = 4.5       # sparse points whose SHAPE encodes a variable, and which
                      # therefore have to be big enough to tell o/s/^/D/v apart


def apply_style():
    """Set every rcParam that controls typeface, size and line weight."""
    plt.rcParams.update({
        # Arial for text and for mathtext, so an axis label and the symbols
        # inside it are the same typeface.
        'font.family': 'sans-serif',
        'font.sans-serif': ['Arial', 'Helvetica', 'Liberation Sans', 'DejaVu Sans'],
        'mathtext.fontset': 'custom',
        'mathtext.rm': 'Arial', 'mathtext.it': 'Arial:italic',
        'mathtext.bf': 'Arial:bold', 'mathtext.default': 'it',
        # one size everywhere
        'font.size': FONT_SIZE, 'axes.labelsize': FONT_SIZE, 'axes.titlesize': FONT_SIZE,
        'xtick.labelsize': FONT_SIZE, 'ytick.labelsize': FONT_SIZE,
        'legend.fontsize': FONT_SIZE, 'legend.title_fontsize': FONT_SIZE,
        # opaque white legend boxes with a 0.5 pt black frame
        'legend.frameon': True, 'legend.facecolor': 'white', 'legend.edgecolor': 'black',
        'legend.framealpha': 1.0, 'legend.fancybox': False,
        'legend.handletextpad': 0.4, 'legend.labelspacing': 0.25,
        'legend.borderpad': 0.35, 'legend.columnspacing': 1.0,
        'patch.linewidth': LW_THIN,
        # lines and markers
        'lines.linewidth': LW_DATA, 'lines.markersize': MS,
        'lines.markeredgewidth': MEW,
        'axes.linewidth': LW_THIN,
        'xtick.major.width': LW_THIN, 'ytick.major.width': LW_THIN,
        'xtick.minor.width': LW_THIN, 'ytick.minor.width': LW_THIN,
        'grid.linewidth': LW_THIN,
        # TrueType, never Type 3 -- ACS production rejects Type 3 subsets
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.facecolor': 'white', 'savefig.transparent': False,
    })


def add_minor_grid(ax):
    """Grid on major and minor ticks, at the 0.5 pt minimum but pale enough to
    stay behind the data.

    Minor grid lines must never land on the exact midpoint between two labelled
    major ticks: a line at 0.015 sitting halfway between 0.01 and 0.02 reads as
    a value the axis never names, and the grid becomes harder to use rather than
    easier.  On LINEAR axes matplotlib's AutoMinorLocator subdivides by 5 whenever
    the major step has mantissa 1, 2.5, 5 or 10, which covers every linear axis
    here, so it is left alone; a major step of 2 would subdivide by 4 and would
    need the same care.  LOG axes are the ones
    that need care, because a hand-written FixedLocator can easily be given
    half-decade positions (1.5, 2.5, 3.5 x 10^n); the call sites here list only
    integer multiples of the decade.
    """
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(which='major', ls='-', lw=LW_THIN, color='0.85')
    ax.grid(which='minor', ls='-', lw=LW_THIN, color='0.93')


def label_panels(axes, dx=-30, dy=2):
    """Bold (a), (b), (c)... just outside the top-left corner of each axes.

    The offset is in POINTS, not axes fractions, so every panel label in the
    paper sits the same distance out no matter how wide its panel is.  As a
    fraction, a wide panel (the 2x2 SI grids) pushed the label off the canvas
    while a narrow one (a 1-of-3 panel) tucked it against the tick labels.
    """
    for ax, letter in zip(axes, 'abcdefgh'):
        ax.annotate(f'({letter})', xy=(0, 1), xycoords='axes fraction',
                    xytext=(dx, dy), textcoords='offset points',
                    fontsize=FONT_SIZE, fontweight='bold', va='bottom', ha='left')


def save_fig(fig, filename):
    """Save as vector PDF at exactly the requested figure size, NEXT TO THE SCRIPT.

    A bare filename is resolved against the directory of the calling script, not
    the current working directory.  Every figure therefore lands beside the code
    that produced it however the script was launched -- `make`, an editor, or a
    shell sitting somewhere else.  Running a script from the wrong directory used
    to scatter stray duplicate PDFs around the tree, which then went stale while
    the real ones were regenerated.  Pass an absolute path to override.

    Deliberately no bbox_inches='tight': the media box has to stay exactly
    SINGLE_COL or DOUBLE_COL wide so the journal reproduces the graphic at
    100 % and the 8 pt lettering prints at 8 pt.  A trimmed box would be
    rescaled at typesetting and the compliance arithmetic would no longer hold.
    Call fig.tight_layout() (or set constrained layout) before this.
    """
    if not os.path.isabs(filename):
        caller = inspect.stack()[1].filename
        filename = os.path.join(os.path.dirname(os.path.abspath(caller)), filename)
    w, h = fig.get_size_inches()
    fig.savefig(filename, format='pdf')
    ok_w = 'single' if abs(w - SINGLE_COL) < 0.01 else \
           'double' if abs(w - DOUBLE_COL) < 0.01 else 'OFF-SPEC'
    print(f'wrote {os.path.basename(filename)}  {w:.3f} x {h:.3f} in ({72 * w:.0f} x {72 * h:.0f} pt)'
          f'  width={ok_w}  depth={"ok" if h <= MAX_DEPTH else "OVER 9.167 in"}')

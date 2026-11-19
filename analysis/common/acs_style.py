"""House figure style, imported by every plotting script.

The numbers come from Appendix 2, "Preparing Graphics", of the ACS author
guide: single-column graphics up to 3.33 in, double-column up to 7.0 in,
lettering at least 4.5 pt in Helvetica or Arial, lines at least 0.5 pt.  The
style was kept for the PNAS submission; see the repository README.  Every
figure is drawn at exactly one of the two column widths and embedded at
natural size, so the 8 pt lettering drawn here is 8 pt in print.
"""

import os

import matplotlib.pyplot as plt
from matplotlib.colors import hsv_to_rgb
from matplotlib.font_manager import FontProperties
from matplotlib.lines import Line2D
from matplotlib.offsetbox import AnchoredOffsetbox, DrawingArea, HPacker, TextArea, VPacker
from matplotlib.textpath import text_to_path
from matplotlib.ticker import FixedFormatter, FixedLocator, NullFormatter

# --- geometry (in) ---
SINGLE_COL = 3.33     # 240 pt
DOUBLE_COL = 7.00     # 504 pt
MAX_DEPTH = 9.167     # 660 pt, including the caption

# --- type ---
FONT_SIZE = 8         # pt, every character in every figure

# --- line weights (pt); nothing below 0.5 ---
LW_THIN = 0.5         # spines, ticks, grids, error bars, guide lines
LW_DATA = 1.0         # data and model curves
MEW = 0.6             # open-marker edge width

# --- markers (pt) ---
MS = 3.5              # the default: I(Q), g2
MS_SPARSE = 4.5       # sparse points whose shape encodes a variable


def apply_style():
    """Set every rcParam that controls typeface, size and line weight."""
    plt.rcParams.update({
        # Arial for text and for mathtext, so a label and its symbols match
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
        # TrueType fonts in the PDF, never Type 3
        'pdf.fonttype': 42, 'ps.fonttype': 42,
        'savefig.facecolor': 'white', 'savefig.transparent': False,
    })


def add_minor_grid(ax):
    """Pale grid on major and minor ticks, behind the data.

    Hand-written log locators in this repository list only integer multiples
    of a decade, so no minor grid line falls halfway between two labelled ticks.
    """
    ax.minorticks_on()
    ax.set_axisbelow(True)
    ax.grid(which='major', ls='-', lw=LW_THIN, color='0.85')
    ax.grid(which='minor', ls='-', lw=LW_THIN, color='0.93')


def q_log_ticks(ax, axis='x'):
    """Label 0.004, 0.01 and 0.03 on a log Q axis.

    The Q range spans about one decade, where matplotlib's own log locator
    labels a single tick.  The same three ticks are used on every Q axis.
    """
    if axis == 'x':
        a = ax.xaxis
    else:
        a = ax.yaxis
    a.set_major_locator(FixedLocator([4e-3, 1e-2, 3e-2]))
    a.set_major_formatter(FixedFormatter(['0.004', '0.01', '0.03']))
    a.set_minor_locator(FixedLocator([3e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                      2e-2, 4e-2, 5e-2]))
    a.set_minor_formatter(NullFormatter())


def label_panels(axes, dx=-30, dy=5, start=0):
    """Bold (a), (b), (c)... just outside the top-left corner of each axes.

    dx and dy are offsets in points, so every label sits the same distance from
    its panel whatever the panel width.  start skips that many letters, for a
    figure whose panels are labelled in two calls.
    """
    letters = 'abcdefghijklmnopqrstuvwxyz'[start:]
    for k in range(len(axes)):
        axes[k].annotate(f'({letters[k]})', xy=(0, 1), xycoords='axes fraction',
                         xytext=(dx, dy), textcoords='offset points',
                         fontsize=FONT_SIZE, fontweight='bold', va='bottom', ha='left')


def text_width(s, prop):
    """Width in points of the string s drawn in the font prop."""
    is_math = s.count('$') >= 2
    return text_to_path.get_text_width_height_descent(s, prop, ismath=is_math)[0]


def pad_to_width(box, width, column_width):
    """box followed by blank space, so the pair is column_width wide."""
    return HPacker(children=[box, DrawingArea(column_width - width, 0)],
                   align='baseline', pad=0, sep=0)


def key_rows(fig, rows, y, x=0.5):
    """One framed key below the panels, with one row per variable.

    rows is a list of (title, handles); each handle is a Line2D carrying the
    marker style and the label.  Each row starts with its title, the titles
    share the first column, and every column is as wide as its widest entry,
    so the markers of all rows line up.  Built from offset boxes rather than
    fig.legend, which would set a title with superscripts on a lower baseline
    than its row.  The box is centred on x with its bottom at y (figure
    fractions).
    """
    fs = FONT_SIZE
    prop = FontProperties(size=fs)
    handle_length = 1.0 * fs
    handle_pad = 0.4 * fs
    column_gap = 1.0 * fs

    # the title of each row, with its width
    titles = []
    for title, handles in rows:
        titles.append((TextArea(title, textprops={'fontsize': fs}), text_width(title, prop)))

    # each entry is a marker in a box the size of a legend handle, then its label
    entries = []
    for title, handles in rows:
        row_entries = []
        for h in handles:
            mark = Line2D([handle_length / 2], [0.35 * fs])
            mark.update_from(h)
            box = DrawingArea(handle_length, 0.7 * fs)
            box.add_artist(mark)
            label = h.get_label()
            packed = HPacker(children=[box, TextArea(label, textprops={'fontsize': fs})],
                             align='baseline', pad=0, sep=handle_pad)
            row_entries.append((packed, handle_length + handle_pad + text_width(label, prop)))
        entries.append(row_entries)

    title_column = 0.0
    for box, width in titles:
        title_column = max(title_column, width)
    n_columns = 0
    for row_entries in entries:
        n_columns = max(n_columns, len(row_entries))
    column_widths = []
    for c in range(n_columns):
        widths = []
        for row_entries in entries:
            if c < len(row_entries):
                widths.append(row_entries[c][1])
        column_widths.append(max(widths))

    lines = []
    for r in range(len(rows)):
        title_box, title_w = titles[r]
        children = [pad_to_width(title_box, title_w, title_column)]
        for c in range(len(entries[r])):
            box, width = entries[r][c]
            children.append(pad_to_width(box, width, column_widths[c]))
        lines.append(HPacker(children=children, align='baseline', pad=0, sep=column_gap))

    key = AnchoredOffsetbox(loc='lower center',
                            child=VPacker(children=lines, align='left', pad=0, sep=0.35 * fs),
                            pad=0.3, borderpad=0, frameon=True, prop={'size': fs},
                            bbox_to_anchor=(x, y), bbox_transform=fig.transFigure)
    key.patch.set_boxstyle('square', pad=0)
    key.patch.set(facecolor='white', edgecolor='black', linewidth=LW_THIN)
    fig.add_artist(key)
    return key


def contrast_on_white(rgb):
    """WCAG contrast ratio of an sRGB colour against a white page."""
    linear = []
    for c in rgb:
        if c <= 0.04045:
            linear.append(c / 12.92)
        else:
            linear.append(((c + 0.055) / 1.055) ** 2.4)
    return 1.05 / (0.2126 * linear[0] + 0.7152 * linear[1] + 0.0722 * linear[2] + 0.05)


def hsv_colors(hues, min_contrast=3.0):
    """Saturated colours at the given hsv hues (0 red, 1/3 green, 2/3 blue),
    each darkened just enough to reach min_contrast against white."""
    colors = []
    for h in hues:
        v = 1.0
        while contrast_on_white(hsv_to_rgb((h, 1.0, v))) < min_contrast:
            v -= 0.01
        colors.append(tuple(hsv_to_rgb((h, 1.0, v))))
    return colors


def save_fig(fig, path):
    """Write fig as a vector PDF at exactly its figure size, and report the size.

    path is where the PDF goes; the scripts pass a path next to themselves.  No
    bbox_inches='tight': the page must stay exactly one column width wide.
    """
    w, h = fig.get_size_inches()
    fig.savefig(path, format='pdf')
    if abs(w - SINGLE_COL) < 0.01:
        width_ok = 'single'
    elif abs(w - DOUBLE_COL) < 0.01:
        width_ok = 'double'
    else:
        width_ok = 'OFF-SPEC'
    if h <= MAX_DEPTH:
        depth_ok = 'ok'
    else:
        depth_ok = 'OVER 9.167 in'
    print(f'wrote {os.path.basename(path)}  {w:.3f} x {h:.3f} in '
          f'({72 * w:.0f} x {72 * h:.0f} pt)  width={width_ok}  depth={depth_ok}')

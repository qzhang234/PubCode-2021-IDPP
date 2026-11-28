"""Figure 2 and SI Figs. S4-S8: the four isothermal experiments of 2021-3.

Writes, next to this script,
  FigureS7-S10_Temperature_*.pdf   SI Figs. S4-S7, one per temperature
                                   (28.00, 28.25, 28.50, 28.75 C)
  Figure3_Isothermal_SAXPCS.pdf    main Fig. 2, the top row of the 28.50 C figure
  FigureS12_Fast_Mode.pdf          SI Fig. S8, the fast mode of every group
and prints the numbers the text quotes with them.  Each per-temperature figure
shows (a) the absolute, background-subtracted I(Q), (b) g2 at the lowest Q with
the fits of Equation 2, (c) the fast fraction f in the five Q bins with D0 on
the right axis, and (d)-(g) g2 at the four other Q bins.

Reads data_2021-3/ only.  The analysis choices (clock, background, contrast,
fit, group selection) are explained in the repository README.
"""

import os
import sys

import h5py
import matplotlib.pyplot as plt
import matplotlib.transforms as mtransforms
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import FuncFormatter, LogLocator

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,  # noqa: E402
                              apply_style, add_minor_grid, hsv_colors, key_rows,
                              label_panels, q_log_ticks, save_fig)
import xpcs_fit                                                           # noqa: E402
from nexus_read import read_start_time, read_g2, read_iq                  # noqa: E402
from abs_xsec import abs_xsec_coef, coef_from_monitors                    # noqa: E402
from series_2021_3 import (DATA, BACKGROUND, BAND, BG_COEF, T_HUES, T_MARKERS,  # noqa: E402
                           Q_BG, I_BG, COEF_BG, abs_iq, acquisitions, band_intensity,
                           group_files, group_midpoint, trace_rows)
from nexus_read import parse_name                                         # noqa: E402

STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')

# header -> (run name in the trace, sample, target C), in order of temperature
RUNS = {'H0402': ('H0402_D8_280C10p', 'D8', 28.00),
        'H0407': ('H0407_D6_282C10p', 'D6', 28.25),
        'H0404': ('H0404_D7_285C10p', 'D7', 28.50),
        'H0413': ('H0413_D9_287C10p', 'D9', 28.75)}
# The 6 C reference of each loading.  D8's is H0403, measured after its
# isothermal experiment: H0402, measured before it, still holds undissolved
# material (README).
REFS = {'H0402': 'H0403_D8_060C10p', 'H0407': 'H0407_D6_060C10p',
        'H0404': 'H0404_D7_060C10p', 'H0413': 'H0413_D9_060C10p'}
QI = [0, 1, 2, 3, 4]                  # the five lowest Q bins, 175 to 83 nm
Q_MARKERS = ['o', 's', '^', 'D', 'v']
ERR_MAX = 0.016                       # largest median g2_err at the lowest Q (1e-4 to 1e-2 s) of a fitted group
CHI2_MAX = 10.0                       # largest reduced chi^2 of a kept fit
CONVERTING = 0.85                     # f at the lowest Q below which arrest is under way
RATE_MIN_F = 0.25                     # below this f at the lowest Q the fast decay fixes no rate
LATE = 0.70                           # f at the lowest Q that separates early and late stages
COLOR_6C = '#1f77b4'                  # the 6 C reference, blue squares in panel (a)
COLOR_20C = '0.45'                    # the first 20 acquisitions of the ramp (~20 C), grey diamonds in (a)

FIG_NAMES = {'H0402': 'FigureS7_Temperature_2800.pdf', 'H0407': 'FigureS8_Temperature_2825.pdf',
             'H0404': 'FigureS9_Temperature_2850.pdf', 'H0413': 'FigureS10_Temperature_2875.pdf'}
MAIN_RUN = 'H0404'                    # 28.50 C, also drawn as Figure 2
MAIN_FIG = 'Figure3_Isothermal_SAXPCS.pdf'
FAST_FIG = 'FigureS12_Fast_Mode.pdf'
MAX_SHOWN = 7                         # more groups than this: (a), (b), (d)-(g) show every second
FIG_T_SIZE = (DOUBLE_COL, 5.7)        # SI figures
FIG_M_SIZE = (DOUBLE_COL, 2.90)       # Figure 2: the top row only
LEGEND_H = 0.072                      # figure fraction kept for the key under the SI panels
G2_YLIM = (0.98, 1.20)

T_COLOR = {}
colors = hsv_colors(T_HUES)
run_keys = list(RUNS)
for k in range(len(run_keys)):
    T_COLOR[run_keys[k]] = colors[k]


def ref_file(h):
    """The averaged 6 C reference of the loading with header h."""
    return os.path.join(DATA, f'Average_{REFS[h]}_att00_Rq0_00001_00020_results.hdf')


BETA = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')
xpcs_fit.CONTRAST = BETA


# --- the background ----------------------------------------------------------
def band_mean(hf):
    """Mean of the stored (Q, phi) intensities inside BAND, the quantity the
    I_band column of the trace holds (no azimuthal average first)."""
    s = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()]).reshape(-1)
    bin_index = hf['/xpcs/qmap/static_index_mapping'][()]
    n_phi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
    q = hf['/xpcs/qmap/static_v_list_dim0'][()][bin_index // n_phi]
    return float(np.nanmean(s[(q >= BAND[0]) & (q < BAND[1])]))


with h5py.File(BACKGROUND, 'r') as hf:
    BAND_BG = band_mean(hf)
print(f'background E0411: coefficient {COEF_BG:.4g}, absolute band mean '
      f'{COEF_BG * BAND_BG:.4f} mm^-1')

# The coefficient that makes each 6 C reference flattest, and the shape that
# BG_COEF leaves.  Each reference is reduced to its median in FIT_NBINS
# logarithmic Q bins and divided by its plateau, so every decade of Q counts
# equally (README, "Background subtraction").
FIT_BAND = (0.0036, 0.0305)           # A^-1; starts at the lowest fitted Q bin
FIT_NBINS = 9
PLATEAU_BAND = (0.012, 0.0305)        # A^-1, where every reference is flat
ALPHA_GRID = np.arange(0.70, 1.15, 0.0025)


def binned_median(q, y):
    """Median of y in each of FIT_NBINS logarithmic Q bins across FIT_BAND."""
    edges = np.exp(np.linspace(np.log(FIT_BAND[0]), np.log(FIT_BAND[1]), FIT_NBINS + 1))
    medians = np.zeros(FIT_NBINS)
    for k in range(FIT_NBINS):
        medians[k] = np.median(y[(q >= edges[k]) & (q < edges[k + 1])])
    return medians


def flatness(q, y):
    """RMS fractional departure of the binned median from the plateau."""
    plateau = np.median(y[(q >= PLATEAU_BAND[0]) & (q < PLATEAU_BAND[1])])
    binned = binned_median(q, y)
    return float(np.sqrt(np.mean((binned / plateau - 1.0) ** 2)))


def shape(q, y):
    """(rise, worst upward step) of the binned median across FIT_BAND.

    rise: lowest bin over highest bin, 1 when flat.  step: largest fractional
    increase from one bin to the next going up in Q, 0 when it never turns up.
    """
    binned = binned_median(q, y)
    step = float(max(0.0, (np.diff(binned) / binned[:-1]).max()))
    return float(binned[0] / binned[-1]), step


def flattest_coef(path):
    """Fraction of the empty cell that makes the 6 C reference in path flattest."""
    with h5py.File(path, 'r') as hf:
        I = abs_xsec_coef(hf) * read_iq(hf)[1]
    cost = []
    for a in ALPHA_GRID:
        cost.append(flatness(Q_BG, I - a * COEF_BG * I_BG))
    return float(ALPHA_GRID[int(np.argmin(cost))])


flattest_text = []
for h in REFS:
    flattest_text.append(f'{h} {flattest_coef(ref_file(h)):.3f}')
print(f'background coefficient: using {BG_COEF:.3f}; flattest would be ' + ', '.join(flattest_text))
for h in REFS:
    with h5py.File(ref_file(h), 'r') as hf:
        I = abs_xsec_coef(hf) * read_iq(hf)[1]
    rise, step = shape(Q_BG, I - BG_COEF * COEF_BG * I_BG)
    print(f'   {REFS[h]}: at {BG_COEF:.2f} the 6 C reference rises {rise:.2f}x toward low Q '
          f'and turns back up by at most {step * 100:.0f}%')


# --- per-acquisition trace: time, temperature, absolute band intensity -------
def load_trace(run):
    """Per acquisition: seconds from t_w = 0, RTD reading (C) and absolute,
    background-subtracted band intensity (mm^-1)."""
    frame, t, T, t0 = acquisitions(run)
    rows = trace_rows(run)
    I = np.zeros(len(rows))
    for k in range(len(rows)):
        coef = coef_from_monitors(float(rows[k]['incident']), float(rows[k]['transmitted']))
        I[k] = coef * float(rows[k]['I_band'])
    I = I - BG_COEF * COEF_BG * BAND_BG
    return t, T, I


# --- fit every group of every isothermal experiment --------------------------
traces = {}
series = {}                           # header -> one dict per fitted group
hold_files = {}                       # header -> every group at the target, fitted or not
for h, (run, sample, target) in RUNS.items():
    traces[h] = load_trace(run)
    frame, t_acq, _, t0 = acquisitions(run)
    rows = []
    hold_files[h] = []
    for fp in group_files(run):
        with h5py.File(fp, 'r') as hf:
            start = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if start < 0:
            continue                  # a group of the ramp
        hold_files[h].append(fp)
        if np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue                  # too noisy to fit
        r = xpcs_fit.fit_seeded(tau, g2, g2e, QI, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        f = np.zeros(len(QI))
        fe = np.zeros(len(QI))
        for k in range(len(QI)):
            f[k] = r['per_q'][QI[k]]['f']
            fe[k] = r['per_q'][QI[k]]['f_err']
        rows.append({'tw': group_midpoint(fp, frame, t_acq), 'fit': r, 'tau': tau, 'g2': g2,
                     'g2e': g2e, 'q': q, 'fp': fp, 'f': f, 'fe': fe})
    series[h] = rows


# --- the numbers the text rests on -------------------------------------------
print(f'\n{"run":19s} {"T":>5s} {"groups":>6s} {"chi2":>9s} {"f0 first->last":>14s} '
      f'{"ordered":>7s} {"t(f0=0.5), group midpoints":>27s}')
HALF = {}
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    f0 = np.array([r['f'][0] for r in rows])
    tw = np.array([r['tw'] for r in rows])
    # groups in which arrest is under way, and how many of them have f rising with Q
    n_arresting = 0
    n_ordered = 0
    for r in rows:
        if r['f'][0] < CONVERTING:
            n_arresting += 1
            if np.all(np.diff(r['f']) > 0):
                n_ordered += 1
    below = np.nonzero(f0 < 0.5)[0]
    if len(below) > 0 and below[0] > 0:
        j = below[0]
        HALF[target] = np.interp(0.5, [f0[j], f0[j - 1]], [tw[j], tw[j - 1]])
        half_text = f'{HALF[target]:8.0f} s'
    else:
        half_text = 'never'
    chi = [r['fit']['red_chi2'] for r in rows]
    print(f'{run} {target:5.2f} {len(rows):6d} {min(chi):4.2f}-{max(chi):4.2f} '
          f'{f0[0]:5.2f}->{f0[-1]:5.2f}   {n_ordered:3d}/{n_arresting:<3d} {half_text:>27s}')
print(f'time to f0 = 0.5 shortens {HALF[28.25] / HALF[28.75]:.1f}-fold from 28.25 to 28.75 C')
span_text = []
for h in RUNS:
    span_text.append(f'{RUNS[h][2]:.2f} C {series[h][0]["tw"]:.0f}-{series[h][-1]["tw"]:.0f} s')
print('first and last t_w (group midpoints): ' + '; '.join(span_text))

print(f'\nD0 (1e7 A^2/s, the 1/e rate) in the early (f0 >= {LATE}) and late (f0 < {LATE}) stages;'
      f'\nthe late range counts only groups with f0 >= {RATE_MIN_F}, where the fast decay fixes a rate:')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    dm = np.array([r['fit']['D0'] for r in rows]) / 1e7
    pf = np.array([r['fit']['p1'] for r in rows])
    f0 = np.array([r['f'][0] for r in rows])
    early = dm[f0 >= LATE]
    late = (f0 < LATE) & (f0 >= RATE_MIN_F)
    if len(early) > 1:
        txt = f'  {run}: early {early[0]:.2f} -> {early[-1]:.2f} ({early[0] / early[-1]:.1f}x)'
    else:
        txt = f'  {run}: early {early[0]:.2f} (one group)'
    if late.sum() > 1:
        txt += (f'; late {dm[late].min():.2f}-{dm[late].max():.2f} ({dm[late].max() / dm[late].min():.1f}x)'
                f' while f0 falls {f0[late].max():.2f} -> {f0[late].min():.2f}')
    unresolved = f0 < RATE_MIN_F
    if unresolved.any():
        txt += f'; unresolved after that: f0 {f0[unresolved].max():.2f} -> {f0[unresolved].min():.2f}'
    print(txt + f'; p_fast {pf.max():.2f} -> {pf.min():.2f}')

print('\nthe same structure, mobile or arresting: I(Q) at the first converting group over')
print('I(Q) of the last 28.00 C group, which is still mobile')
q_ref, I_ref = abs_iq(series['H0402'][-1]['fp'])
for h in ['H0407', 'H0404', 'H0413']:
    first = None
    for r in series[h]:
        if r['f'][0] < CONVERTING:
            first = r
            break
    q, I = abs_iq(first['fp'])
    ratio_text = []
    for x in (0.004, 0.006, 0.008, 0.012, 0.02, 0.028):
        ratio_text.append(f'{np.interp(x, q, I) / np.interp(x, q_ref, I_ref):.2f}')
    print(f'  {RUNS[h][0]} t_w {first["tw"]:.0f} s, f0 {first["f"][0]:.2f}: ' + ' '.join(ratio_text))

print('\nstructure against arrest (absolute band intensity of the group I(Q), mm^-1,'
      ' as Figure 3c draws it):')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    arresting = []
    for r in rows:
        if r['f'][0] < CONVERTING:
            arresting.append(r)
    if len(arresting) == 0:
        I_end = band_intensity(rows[-1]['fp'])
        I_0 = band_intensity(hold_files[h][0])
        f0_min = np.min([r['f'][0] for r in rows])
        print(f'  {run}: no arrest; band {I_0:.2f} over the first group -> {I_end:.2f} at {rows[-1]["tw"]:.0f} s '
              f'({I_end / I_0:.0f}x), f0 >= {f0_min:.2f} throughout')
        continue
    a = arresting[0]
    b = arresting[-1]
    Ia = band_intensity(a['fp'])
    Ib = band_intensity(b['fp'])
    qa, ia = abs_iq(a['fp'])
    qb, ib = abs_iq(b['fp'])
    ratio_text = []
    for x in (0.004, 0.008, 0.02):
        ratio_text.append(f'{x:.3f}:{np.interp(x, qa, ib) / np.interp(x, qa, ia):.2f}')
    print(f'  {run}: {a["tw"]:.0f} -> {b["tw"]:.0f} s: 1-f0 {1 - a["f"][0]:.2f} -> {1 - b["f"][0]:.2f}, '
          f'band {Ia:.2f} -> {Ib:.2f} ({Ib / Ia:.2f}x); I(Q) ratio last/first at Q ' + '  '.join(ratio_text))

# How far the structure has grown when the dynamics start to arrest.  The ~20 C
# level is the first 20 acquisitions of the ramp, the reference
# arrest_vs_temperature.py quotes the half-arrest intensity against.
print('\nlow-Q band intensity of every group on target (mm^-1), with f at the lowest Q '
      'where fitted,\nagainst the first 20 acquisitions of the ramp (~20 C):')
for h, (run, sample, target) in RUNS.items():
    frame, t_acq, _, _ = acquisitions(run)
    base20 = band_intensity(group_files(run)[0])
    fitted_f0 = {}
    for r in series[h]:
        fitted_f0[r['fp']] = r['f'][0]
    cells = []
    for fp in hold_files[h]:
        cell = f'{group_midpoint(fp, frame, t_acq):.0f} s {band_intensity(fp):.2f}'
        if fp in fitted_f0:
            cell += f' (f0 {fitted_f0[fp]:.2f})'
        cells.append(cell)
    print(f'  {target:.2f} C, ~20 C level {base20:.4f}: ' + '; '.join(cells))
    cells = []
    for fp in hold_files[h]:
        cells.append(f'{group_midpoint(fp, frame, t_acq):.0f} s {band_intensity(fp) / base20:.0f}x')
    print('     ' + ', '.join(cells) + ' the ~20 C level')

print('\nramp: band intensity relative to the 20-22 C level of each loading')
for h, (run, sample, target) in RUNS.items():
    t, T, I = traces[h]
    base = np.median(I[T < 22.0])
    cells = []
    for x in (26.0, 27.0, 27.5, 27.75, 28.0):
        if x < target:
            near = (T >= x - 0.13) & (T < x + 0.13) & (t < 0)
            cells.append(f'{x:.2f}C:{np.median(I[near]) / base:.1f}')
    print(f'  {run}: 20-22 C {base:.4f} mm^-1; ' + '  '.join(cells))

print('\n6 C references, absolute and background-subtracted (mm^-1):')
for h in RUNS:
    q, I = abs_iq(ref_file(h))
    flat = (q > 0.006) & (q < 0.03)
    print(f'  {REFS[h]}: mean over 0.006-0.03 {np.mean(I[flat]):.4f} +/- {np.std(I[flat]):.4f}')

# How far the empty cell stands above the 6 C references, and how far the
# fitted groups stand above the subtracted background: the two ratios of
# SI section 6.2.
print('\nempty cell against the 6 C reference it has to be subtracted from:')
bg_all = COEF_BG * I_BG               # all of the empty cell, absolute
for q_target in (0.0036, 0.0040, 0.0050):
    j = int(np.argmin(np.abs(Q_BG - q_target)))
    ratios = []
    for h in RUNS:
        I_6C = abs_iq(ref_file(h))[1]
        ratios.append(bg_all[j] / I_6C[j])
    print(f'  Q = {Q_BG[j]:.5f}: empty cell {bg_all[j]:.4g} mm^-1, '
          f'{min(ratios):.1f}-{max(ratios):.1f}x the 6 C signal '
          f'({BG_COEF:.2f} x it: {min(ratios)*BG_COEF:.1f}-{max(ratios)*BG_COEF:.1f}x)')

print('\nfitted groups against the subtracted background, over 0.004-0.008 A^-1:')
in_band = (Q_BG >= 0.004) & (Q_BG <= 0.008)
bg_used = BG_COEF * COEF_BG * I_BG    # the part of the empty cell subtracted
all_points = []
all_means = []
for h, (run, sample, target) in RUNS.items():
    points = []
    means = []
    for r in series[h]:
        ratio = abs_iq(r['fp'])[1][in_band] / bg_used[in_band]
        points.append(ratio)
        means.append(ratio.mean())
    all_points += points
    all_means += means
    pointwise = np.concatenate(points)
    print(f'  {run}: {len(means):2d} groups, pointwise {pointwise.min():5.0f}-{pointwise.max():5.0f}x, '
          f'group means {min(means):5.0f}-{max(means):5.0f}x')
pointwise = np.concatenate(all_points)
print(f'  all {len(all_means)} fitted groups: pointwise {pointwise.min():.0f}-{pointwise.max():.0f}x, '
      f'group means {min(all_means):.0f}-{max(all_means):.0f}x')
print('  (SI section 6.2 quotes the group-mean range; the two differ because one is a '
      'per-point\n   extremum and the other a per-group average)')


# --- Figure 2 and Figures S4-S7 ----------------------------------------------
apply_style()


def every_other_decade(x, position):
    """Tick label 10^n for even n, blank for odd n."""
    power = int(round(np.log10(x)))
    if power % 2 == 0:
        return rf'$10^{{{power}}}$'
    return ''


def time_color(tw, t_lo, t_hi):
    """Colour of a waiting time: plasma, from t_lo to t_hi."""
    if t_hi == t_lo:
        x = 0.5
    else:
        x = (tw - t_lo) / (t_hi - t_lo)
    return plt.cm.plasma(0.10 + 0.78 * x)


for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    if len(rows) <= MAX_SHOWN:
        shown = rows
    else:
        shown = rows[::2]
    if shown[-1] is not rows[-1]:
        shown = shown + [rows[-1]]

    # Every group at the target is drawn in (a), fitted or not, so the colour
    # scale spans all of them; (c) keeps the range of the fitted groups.
    frame, t_acq, T_acq, _ = acquisitions(run)
    hold_tw = []
    for fp in hold_files[h]:
        hold_tw.append(group_midpoint(fp, frame, t_acq))
    t_lo = min(hold_tw)
    t_hi = max(hold_tw)
    fitted_files = []
    for r in rows:
        fitted_files.append(r['fp'])
    shown_files = []
    for r in shown:
        shown_files.append(r['fp'])

    # The first 20 acquisitions of the ramp: the ~20 C level that the text
    # quotes the low-Q growth against.
    ramp_first = group_files(run)[0]
    _, a20, b20 = parse_name(ramp_first)
    in20 = (frame >= a20) & (frame <= b20)
    print(f'{target:.2f} C: ~20 C reference = acquisitions {a20}-{b20}, RTD {T_acq[in20].min():.2f}-'
          f'{T_acq[in20].max():.2f} C; {len(hold_files[h])} groups at the target, {len(rows)} fitted')
    print(f'   acquisitions start {-t_acq[0] / 60:.1f} min before t_w = 0, at RTD {T_acq[0]:.2f} C, '
          f'{int(np.sum(t_acq < 0))} of them on the ramp')

    # The 28.50 C loading is drawn twice: the seven-panel SI figure and the
    # three-panel Figure 2.
    if h == MAIN_RUN:
        versions = [False, True]
    else:
        versions = [False]
    for main_text in versions:
        if main_text:
            fig = plt.figure(figsize=FIG_M_SIZE)
            top = GridSpec(1, 3, figure=fig, left=0.075, right=0.925, top=0.90,
                           bottom=0.325, wspace=0.42)
        else:
            fig = plt.figure(figsize=FIG_T_SIZE)
            top = GridSpec(1, 3, figure=fig, left=0.075, right=0.925, top=0.955,
                           bottom=0.595, wspace=0.42)
        ax1 = fig.add_subplot(top[0, 0])
        ax2 = fig.add_subplot(top[0, 1])
        ax3 = fig.add_subplot(top[0, 2])
        axq = []
        if not main_text:
            bottom = GridSpec(1, 4, figure=fig, left=0.075, right=0.925, top=0.470,
                              bottom=0.088 + LEGEND_H, wspace=0.08)
            for k in range(4):
                axq.append(fig.add_subplot(bottom[0, k]))
            for a in axq[1:]:
                a.sharey(axq[0])
        label_panels([ax1, ax2, ax3])
        if axq:
            label_panels(axq, dx=-14, start=3)   # no y tick labels to clear there

        # (a) SAXS, from the lowest fitted Q bin outward; below it the points
        # sit on the edge of the beamstop shadow (README).  The 6 C reference,
        # the first 20 acquisitions of the ramp (~20 C) and every group at the target:
        # circles where g2 was fitted, triangles where it was too noisy to fit
        # (those appear in (a) only).  The key lists the triangles and the
        # fitted groups drawn in (b) and (d)-(g); the profiles between them take
        # their colour from the same t_w scale.
        key = []
        q6, I6 = abs_iq(ref_file(h))
        m = (I6 > 0) & (q6 >= FIT_BAND[0])
        ax1.plot(q6[m], I6[m], 's', color=COLOR_6C, ms=MS, mfc='none', mew=MEW)
        lo = I6[m].min()
        hi = I6[m].max()
        q20, I20 = abs_iq(ramp_first)
        m = (I20 > 0) & (q20 >= FIT_BAND[0])
        ax1.plot(q20[m], I20[m], 'D', color=COLOR_20C, ms=MS - 0.5, mfc='none', mew=MEW)
        lo = min(lo, I20[m].min())
        hi = max(hi, I20[m].max())
        for k in range(len(hold_files[h])):
            fp = hold_files[h][k]
            c = time_color(hold_tw[k], t_lo, t_hi)
            if fp in fitted_files:
                mk = 'o'
            else:
                mk = '^'
            q, I = abs_iq(fp)
            m = (I > 0) & (q >= FIT_BAND[0])
            ax1.plot(q[m], I[m], mk, color=c, ms=MS, mfc='none', mew=MEW)
            lo = min(lo, I[m].min())
            hi = max(hi, I[m].max())
            if fp in shown_files or fp not in fitted_files:
                key.append(Line2D([], [], marker=mk, ls='none', color=c, ms=MS, mfc='none',
                                  mew=MEW, label=f'{hold_tw[k]:.0f}'))
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(3.3e-3, 3.6e-2)
        q_log_ticks(ax1)
        ax1.set_ylim(lo * 0.3, hi * 6.0)
        # The low-Q band of Figure 3c.  Its edges are dashed lines at zorder 1,
        # above the minor grid lines that fall on the same Q values.
        ax1.axvspan(BAND[0], BAND[1], color='0.92', lw=0, zorder=0)
        for edge in BAND:
            ax1.axvline(edge, color='0.45', ls='--', lw=LW_THIN, zorder=1.0)
        ax1.text(np.sqrt(BAND[0] * BAND[1]), 0.97, 'low-$Q$\nband', ha='center', va='top',
                 color='0.35', linespacing=1.0,
                 transform=mtransforms.blended_transform_factory(ax1.transData, ax1.transAxes))
        plateau = np.median(I6[(q6 > 0.006) & (q6 < 0.03)])
        ax1.text(0.012, plateau * 0.6, '6 °C ref', ha='center', va='top', color=COLOR_6C)
        ax1.text(0.012, np.interp(0.012, q20, I20) * 1.5, '20 °C', ha='center', va='bottom',
                 color=COLOR_20C)
        ax1.set_xlabel(r'$Q$ ($\AA^{-1}$)')
        ax1.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
        add_minor_grid(ax1)

        # (b) and (d)-(g): g2 with the fits
        g2_panels = [(ax2, 0)]
        for k in range(len(axq)):
            g2_panels.append((axq[k], k + 1))
        for ax, qi in g2_panels:
            for r in shown:
                c = time_color(r['tw'], t_lo, t_hi)
                m = r['tau'] > 0
                ax.errorbar(r['tau'][m], r['g2'][m, qi], yerr=r['g2e'][m, qi], fmt='o', color=c,
                            mfc='none', ms=MS - 0.5, mew=LW_THIN, capsize=1.0, elinewidth=LW_THIN,
                            capthick=LW_THIN, alpha=0.75, zorder=2)
                fit = r['fit']
                pq = fit['per_q'][qi]
                tt = np.logspace(np.log10(r['tau'][m].min()), np.log10(r['tau'][m].max()), 200)
                ax.plot(tt, xpcs_fit.double_exp(tt, pq['tau_fast'], pq['f'], pq['tau_slow'],
                                                fit['p1'], fit['p2']), color=c, lw=LW_DATA, zorder=3)
            # 1 + beta, the contrast every fit uses (the plateau of the static standard)
            ax.axhline(1.0 + BETA, color='0.35', ls='--', lw=LW_DATA * 1.7, zorder=1)
            ax.set_xscale('log')
            ax.set_xlim(9e-6, 5.0)
            ax.set_ylim(*G2_YLIM)
            ax.set_yticks([1.00, 1.05, 1.10, 1.15])
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
            ax.xaxis.set_major_formatter(FuncFormatter(every_other_decade))
            ax.set_xlabel(r'Delay Time, $\Delta t$ (s)')
            add_minor_grid(ax)
            ax.text(0.96, 0.95, rf'$Q = {rows[0]["q"][qi]:.5f}\ \AA^{{-1}}$', transform=ax.transAxes,
                    ha='right', va='top')
        ax2.annotate(r'$1+\beta$', xy=(0.985, 1.0 + BETA),
                     xycoords=ax2.get_yaxis_transform(), xytext=(0, 4),
                     textcoords='offset points', ha='right', va='bottom', color='0.35')
        ax2.set_ylabel(r'$g_2$')
        if axq:
            axq[0].set_ylabel(r'$g_2$')
            for a in axq[1:]:
                plt.setp(a.get_yticklabels(), visible=False)

        # (c) f against t_w: colour = t_w, marker = Q bin; D0 on the right axis
        for k in range(len(QI)):
            i = QI[k]
            mk = Q_MARKERS[k]
            x = [r['tw'] for r in rows]
            ax3.plot(x, [r['f'][i] for r in rows], '-', color='0.75', lw=LW_THIN, zorder=1)
            for r in rows:
                ax3.errorbar(r['tw'], r['f'][i], yerr=r['fe'][i], marker=mk,
                             color=time_color(r['tw'], t_lo, t_hi),
                             ms=MS_SPARSE - 1, capsize=1.0, mfc='none', mew=MEW,
                             elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
        f_lo = np.inf
        f_hi = -np.inf
        for r in rows:
            f_lo = min(f_lo, np.min(r['f'] - r['fe']))
            f_hi = max(f_hi, np.max(r['f'] + r['fe']))
        span = max(f_hi - f_lo, 0.2)
        ax3.set_ylim(f_lo - 0.08 * span, f_hi + 0.42 * span)
        if span > 0.5:
            step = 0.2
        else:
            step = 0.1
        ticks = []
        for v in np.arange(0.0, 1.0001, step):
            if f_lo - 0.08 * span <= v:
                ticks.append(v)
        ax3.set_yticks(ticks)
        fit_lo = rows[0]['tw']
        fit_hi = rows[-1]['tw']
        pad = 0.12 * max(fit_hi - fit_lo, 600)
        ax3.set_xlim(fit_lo - pad, fit_hi + pad)
        ax3.set_xlabel(r'$t_w$ (s)')
        ax3.set_ylabel('Fast Fraction, $f$')
        add_minor_grid(ax3)
        q_handles = []
        for k in range(len(QI)):
            q_handles.append(Line2D([], [], marker=Q_MARKERS[k], ls='none', mfc='none', mec='k',
                                    mew=MEW, ms=MS_SPARSE - 1,
                                    label=f'{1e3 * rows[0]["q"][QI[k]]:.2f}'))
        # D0: black where the fast decay fixes a rate, open grey where f at the
        # lowest Q is below RATE_MIN_F
        ax3r = ax3.twinx()
        d0_top = 0.0
        for r in rows:
            y = r['fit']['D0'] / 1e7
            e = r['fit']['D0_err'] / 1e7
            if r['f'][0] >= RATE_MIN_F:
                star_color = 'k'
            else:
                star_color = '0.6'
            ax3r.errorbar(r['tw'], y, yerr=e, marker='*', ls='none', color=star_color,
                          ms=MS_SPARSE + 2.0, capsize=1.0, mfc='none', mew=MEW, elinewidth=LW_THIN,
                          capthick=LW_THIN, zorder=4)
            d0_top = max(d0_top, y + e)
        ax3r.set_ylim(0, d0_top * 1.30)
        ax3r.set_ylabel(r'$D_0$ ($10^{7}\ \AA^{2}\,$s$^{-1}$)')

        # one key box under the panels: t_w above Q
        key_rows(fig, [(r'$t_w$ (s)', key), (r'$Q$ ($10^{-3}\ \AA^{-1}$)', q_handles)], y=0.004)
        if main_text:
            save_fig(fig, os.path.join(HERE, MAIN_FIG))
        else:
            save_fig(fig, os.path.join(HERE, FIG_NAMES[h]))


# --- Figure S8: f and the fast mode of every group, against the waiting time --
# Plotted against t_w, which only increases, so each line runs forward in time.
def draw_against_time(ax, rows, y, e, color, marker, open_low_f):
    """One quantity of every group of one temperature against t_w, joined in
    time order; open symbols where f at the lowest Q is below RATE_MIN_F if
    open_low_f is True (the fast decay then fixes no rate)."""
    t = np.array([r['tw'] for r in rows])
    ax.plot(t, y, '-', color=color, lw=LW_THIN, alpha=0.6, zorder=1)
    for k in range(len(rows)):
        if open_low_f and rows[k]['f'][0] < RATE_MIN_F:
            face = 'none'
        else:
            face = color
        ax.errorbar(t[k], y[k], yerr=e[k], marker=marker, ls='none', color=color, ms=MS_SPARSE - 1,
                    mfc=face, mew=MEW, elinewidth=LW_THIN, capsize=1.0, capthick=LW_THIN, zorder=2)


fig = plt.figure(figsize=(DOUBLE_COL, 2.75))
gs = GridSpec(1, 3, figure=fig, left=0.068, right=0.99, top=0.92, bottom=0.30, wspace=0.36)
axf = fig.add_subplot(gs[0, 0])
axd = fig.add_subplot(gs[0, 1])
axp = fig.add_subplot(gs[0, 2])
for k in range(len(run_keys)):
    h = run_keys[k]
    rows = series[h]
    f0 = np.array([r['f'][0] for r in rows])
    f0_err = np.array([r['fe'][0] for r in rows])
    d0 = np.array([r['fit']['D0'] for r in rows]) / 1e7
    d0_err = np.array([r['fit']['D0_err'] for r in rows]) / 1e7
    pf = np.array([r['fit']['p1'] for r in rows])
    pf_err = np.array([r['fit']['p1_err'] for r in rows])
    draw_against_time(axf, rows, f0, f0_err, T_COLOR[h], T_MARKERS[k], False)
    draw_against_time(axd, rows, d0, d0_err, T_COLOR[h], T_MARKERS[k], True)
    draw_against_time(axp, rows, pf, pf_err, T_COLOR[h], T_MARKERS[k], True)
# Reference lines sit at zorder 1, above the grid (zorder 0.5), which hid them.
axf.axhline(RATE_MIN_F, color='0.35', ls='--', lw=LW_THIN, zorder=1)
axf.text(5000, RATE_MIN_F - 0.02, f'$f$ = {RATE_MIN_F}', color='0.35', ha='center', va='top')
axf.text(0.97, 0.97, f'$Q$ = {series["H0404"][0]["q"][0]:.5f} ' + r'$\AA^{-1}$',
         transform=axf.transAxes, ha='right', va='top')
axf.set_ylim(-0.08, 1.15)
axf.set_ylabel(r'Fast Fraction, $f$')
axd.set_yscale('log')
axd.set_ylabel(r'$D_0$ ($10^{7}\ \AA^{2}\,$s$^{-1}$)')
axp.axhline(xpcs_fit.P_EXP_LO[0], color='0.35', ls='--', lw=LW_THIN, zorder=1)
axp.text(8700, xpcs_fit.P_EXP_LO[0] - 0.015, 'lower bound of the fit', color='0.35', ha='right',
         va='top')
axp.set_ylim(0.0, 1.1)
axp.set_ylabel(r'$p_{\mathrm{fast}}$')
for ax in (axf, axd, axp):
    ax.set_xlim(-400, 8900)
    ax.set_xlabel(r'$t_w$ (s)')
    add_minor_grid(ax)
label_panels([axf, axd, axp])
handles = [Line2D([], [], ls='none', marker='', label=r'$T$ (°C)')]
for k in range(len(run_keys)):
    h = run_keys[k]
    handles.append(Line2D([], [], ls='none', marker=T_MARKERS[k], color=T_COLOR[h], mfc=T_COLOR[h],
                          mew=MEW, ms=MS_SPARSE - 1, label=f'{RUNS[h][2]:.2f}'))
handles.append(Line2D([], [], ls='none', marker='o', color='0.3', mfc='none', mew=MEW,
                      ms=MS_SPARSE - 1, label=f'open: $f$ < {RATE_MIN_F}, no resolved rate'))
legend = fig.legend(handles=handles, loc='lower center', bbox_to_anchor=(0.5, 0.005),
                    ncol=len(handles), handlelength=1.0, handletextpad=0.35, columnspacing=1.0,
                    borderpad=0.3)
legend.get_frame().set_linewidth(LW_THIN)
save_fig(fig, os.path.join(HERE, FAST_FIG))

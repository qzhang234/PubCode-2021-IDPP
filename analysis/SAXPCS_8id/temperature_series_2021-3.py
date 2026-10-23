"""The Figure 3 analysis applied to the 2021-3 temperature series.

Reads only data_2021-3/, written by reduce_babnigg202110.py.  Writes a CANDIDATE
figure, Candidate_TemperatureSeries_2021-3.pdf, which is not yet part of the
paper.

THE DATA.  Four ELP loadings in one cell block, each taken from 20 C to a final
temperature on the QNW zone-3 controller at 0.10 C/min and then held there:
D8 at 28.00, D6 at 28.25, D7 at 28.50 and D9 at 28.75 C.  t_w = 0 is the first
acquisition with zone 3 at its target, the analogue of the RTD reaching 30 C in
Figure 3.  Temperatures are controller readings, not the calibrated RTD of the
2022-1 series, so they order the four runs but do not place them on the 2022-1
temperature scale.

WHAT IS THE SAME AS FIGURE 3.  The reduction (100-acquisition blocks through the
two outlier cuts of average_ranges.py), the model (Equation 2, the fast mode
tied to 1/(D0 Q^2), beta fixed), the five lowest Q bins, the error-weighted
least squares and covariance, the chi^2 guard of 10, and the absolute scale:
every I(Q) is converted with abs_xsec.py from its own ion-chamber readings (the
2021-3 beamtime used the same chambers) and the background is subtracted as
coef_sample * I_sample - coef_bg * I_bg, exactly as Figure 3 subtracts the buffer.

WHAT DIFFERS, AND WHY.
* The background is the empty cell E0411, not a buffer, because 2021-3 has no
  buffer measurement.  Of the six empty-cell and air measurements it is the only
  one with windows in the beam (transmission 0.94; G0409, G0410 and the air shot
  transmit 1.00), and its low-Q excess over the windowless ones has the shape of
  the low-Q excess of the dissolved samples, at 0.7-1.0 of its amplitude.  With
  it the dissolved state is flat from 0.006 to 0.03 A^-1.  It slightly
  over-subtracts below about 0.0055 A^-1.  E0008 and E0304 carry no transmitted-
  beam monitor and cannot be put on the absolute scale at all.
* One detector pixel, (y, x) = (63, 166), is stuck at 1.5 photons per frame in
  every acquisition, sample and empty alike, and the CrossRm qmap does not mask
  it.  It lies in the static bin at Q = 0.0214 A^-1, phi sector 7, far from the
  five dynamic bins, so no correlation function is affected; its (Q, phi) entry
  is left out of the azimuthal average of I(Q) here.
* beta is measured for THIS beamtime, on the 12IDB static standard B0005, in
  the same five bins and by the same procedure as contrast_calibration.py.  The
  samples' own g2 - 1 below 60 us averages 0.132 over the blocks where arrest is
  under way, against beta = 0.133.
* p_fast is fitted in each block, not shared across the series.  The paper
  shares it because the 2022-1 chi^2 profile is flat between 0.55 and 0.70.
  These runs also capture the mobile stage before arrest, and across it the
  free p_fast moves from about 0.9 to about 0.3; one value cannot describe both.
  Each block therefore gets its own p_fast, started from the minimum of a p_fast
  profile because a single start can land in a degenerate solution where the
  fast term mimics the slow one.
* The fast-mode rate is reported as D_med = 1/(tau_med Q^2), from the MEDIAN
  Kohlrausch time tau_med = tau (ln 2)^(1/p), not as the paper's mean-time D_eff.
  With p_fast free in each block the mean time Gamma(1/p) tau / p is dominated
  by the tail and swings with p: on the Figure 3 groups themselves, freeing
  p_fast widens the D_eff span from 2.9x to 5.0x while D_med spans 1.7x.
* A block is fitted when its correlation function is no noisier than the
  noisiest group Figure 3 itself fits: median g2_err at the lowest Q over
  1e-4 to 1e-2 s of at most 0.016.

Panels: (a) I(Q) through the 28.50 C hold, with the 6 C reference; (b) the mean
I(Q) over 0.004-0.008 A^-1 against zone-3 temperature on the ramp; (c) the same
through the hold; (d) g2 at the lowest Q through the 28.50 C hold, with fits;
(e)-(h) fast fraction f against t_w in the five Q bins, one panel per
temperature, with D_med on the right axis, as Figure 3c carries D_eff.
"""

import csv
import glob
import os
import sys
from datetime import datetime

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
import xpcs_fit                                                  # noqa: E402
from nexus_read import parse_name, read_start_time, read_g2      # noqa: E402
from abs_xsec import (abs_xsec_coef, CAL_A, CAL_B, T_EXP,         # noqa: E402
                      AIR_TRANSMISSION, NUM_FRAMES, SAMPLE_THICKNESS,
                      DELTA_OMEGA)

DATA = os.path.join(HERE, 'data_2021-3')
TRACE = os.path.join(DATA, 'trace_2021-3.csv')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
BACKGROUND = os.path.join(DATA, 'Average_E0411_Empty_3_006C_att00_Rq0_00001_00050_results.hdf')
STUCK_PIXEL = (63, 166)               # (y, x); see the module docstring
BAND = (0.004, 0.008)                 # A^-1, the band of Figure S4b and of the trace

# header -> (run name in the trace, sample, target C), in order of temperature
RUNS = {'H0402': ('H0402_D8_280C10p', 'D8', 28.00),
        'H0407': ('H0407_D6_282C10p', 'D6', 28.25),
        'H0404': ('H0404_D7_285C10p', 'D7', 28.50),
        'H0413': ('H0413_D9_287C10p', 'D9', 28.75)}
REF_6C = os.path.join(DATA, 'Average_H0404_D7_060C10p_att00_Rq0_00001_00020_results.hdf')
QI = [0, 1, 2, 3, 4]                  # the five Q bins, as in Figure 3
Q_MARKERS = ['o', 's', '^', 'D', 'v']
T_MARKERS = ['o', 's', '^', 'D']
T_CMAP, T_LO, T_HI = plt.cm.viridis, 0.05, 0.85
AT_TARGET = 0.01                      # C
ERR_MAX = 0.016                       # the noisiest group Figure 3 fits
CHI2_MAX = 10.0                       # the Figure 3 guard
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
SHOW_RUN = 'H0404'                    # the run shown in panels (a) and (d)
COLOR_6C = '#1f77b4'                  # the 6 C reference, as in Figure 3a
FIG_SIZE = (DOUBLE_COL, 5.0)
D_TOP = 15.0                          # right-axis limit for D_med, 1e7 A^2/s

T_COLOR = {h: T_CMAP(c) for h, c in zip(RUNS, np.linspace(T_LO, T_HI, len(RUNS)))}
INC = '/entry/instrument/incident_beam/incident_beam_intensity'
TRN = '/entry/instrument/incident_beam/transmitted_beam_intensity'


# --- beta for this beamtime, exactly as contrast_calibration.py measures it ---
def flat_level(tau, g, e):
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return np.sum(w * g[ok]) / np.sum(w)


with h5py.File(STANDARD, 'r') as hf:
    tau, g2, g2e, _ = read_g2(hf)
betas = np.array([flat_level(tau, g2[:, i], g2e[:, i]) - 1.0 for i in QI])
BETA = betas.mean()
print(f'beta (2021-3, 12IDB standard, bins 0-4) = {BETA:.5f} +/- '
      f'{betas.std(ddof=1) / np.sqrt(len(betas)):.5f}; per bin '
      + ', '.join(f'{b:.4f}' for b in betas))
xpcs_fit.CONTRAST = BETA              # double_exp reads it at call time


# --- absolute scale ---
def coef_from(up_ic, dn_ic):
    """abs_xsec_coef for monitor readings given as numbers rather than an open file."""
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    transmission = dn_ic / up_ic / AIR_TRANSMISSION
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


def read_iq(hf):
    """Azimuthally averaged I(Q) with the stuck pixel's (Q, phi) entry left out.

    Otherwise the same as nexus_read.read_saxs_iq: one value per (Q, phi) bin,
    flattened, phi running fastest, averaged over phi at each Q.
    """
    intensity = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()], dtype=float).reshape(-1)
    idx = hf['/xpcs/qmap/static_index_mapping'][()]
    q_list = hf['/xpcs/qmap/static_v_list_dim0'][()]
    n_phi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
    stuck_roi = int(hf['/xpcs/qmap/static_roi_map'][STUCK_PIXEL])
    if stuck_roi > 0:
        intensity[stuck_roi - 1] = np.nan          # roi ids count entries from 1
    q_idx = idx // n_phi
    uq = np.unique(q_idx)
    return q_list[uq], np.array([np.nanmean(intensity[q_idx == i]) for i in uq])


def band_mean(hf):
    """Mean of the stored (Q, phi) entries in BAND, as the trace records it."""
    s = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()]).reshape(-1)
    idx = hf['/xpcs/qmap/static_index_mapping'][()]
    q = hf['/xpcs/qmap/static_v_list_dim0'][()][idx // hf['/xpcs/qmap/static_v_list_dim1'].shape[0]]
    return float(np.nanmean(s[(q >= BAND[0]) & (q < BAND[1])]))


with h5py.File(BACKGROUND, 'r') as hf:
    Q_BG, I_BG = read_iq(hf)
    COEF_BG = abs_xsec_coef(hf)
    BAND_BG = band_mean(hf)
    check = coef_from(float(np.asarray(hf[INC][()]).ravel()[0]),
                      float(np.asarray(hf[TRN][()]).ravel()[0]))
assert abs(check / COEF_BG - 1) < 1e-12, 'coef_from disagrees with abs_xsec_coef'
print(f'background E0411: coefficient {COEF_BG:.4g}, absolute band mean '
      f'{COEF_BG * BAND_BG:.4f} mm^-1')


def abs_iq(fp):
    """Absolute, background-subtracted I(Q) of one averaged file, in mm^-1."""
    with h5py.File(fp, 'r') as hf:
        q, I = read_iq(hf)
        coef = abs_xsec_coef(hf)
    assert np.allclose(q, Q_BG), 'sample and background on different Q grids'
    return q, coef * I - COEF_BG * I_BG


# --- per-acquisition trace: absolute band intensity, temperature, time ---
def load_trace(run):
    rows = [r for r in csv.DictReader(open(TRACE)) if r['run'] == run]
    rows.sort(key=lambda r: int(r['frame']))
    when = [datetime.strptime(r['start_time'], '%Y-%m-%d %H:%M:%S') for r in rows]
    T = np.array([float(r['qnw3_C']) for r in rows])
    target = float(rows[0]['qnw3_set_C'])
    I = np.array([coef_from(float(r['incident']), float(r['transmitted'])) * float(r['I_band'])
                  for r in rows]) - COEF_BG * BAND_BG
    i0 = int(np.argmax(T >= target - AT_TARGET))
    t = np.array([(w - when[i0]).total_seconds() for w in when])
    return t, T, I, when[i0]


def binned(x, y, edges):
    """Median of y in each bin of x, for bins holding at least five points."""
    xc, yc = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() >= 5:
            xc.append(np.median(x[m]))
            yc.append(np.median(y[m]))
    return np.array(xc), np.array(yc)


# --- Equation 2 in one block, p_fast free, started at its profile minimum ---
def fit_block(tau, g2, g2e, q):
    prof = [(r['red_chi2'], pf, r) for pf in P_GRID
            for r in [xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=pf)]
            if r is not None]
    chi0, pf0, best = min(prof, key=lambda x: x[0])
    xpcs_fit.P_EXP_P0[0] = pf0
    free = xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q)
    xpcs_fit.P_EXP_P0[0] = 0.5
    if free is not None and free['red_chi2'] <= chi0:
        best = free
    return best


def dmed(fit):
    """Median-time diffusion coefficient: tau_med = tau (ln 2)^(1/p), so D_med = D0 / (ln 2)^(1/p)."""
    return fit['D0'] / np.log(2) ** (1.0 / fit['p1'])


traces, series = {}, {}
for h, (run, sample, target) in RUNS.items():
    t, T, I, t0 = load_trace(run)
    traces[h] = (t, T, I)
    rows = []
    for fp in sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')),
                     key=lambda p: parse_name(p)[1]):
        with h5py.File(fp, 'r') as hf:
            tw = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        noisy = np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX
        if tw < 0 or noisy:
            continue
        r = fit_block(tau, g2, g2e, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        rows.append({'tw': tw, 'fit': r, 'tau': tau, 'g2': g2, 'g2e': g2e, 'q': q, 'fp': fp,
                     'f': np.array([r['per_q'][i]['f'] for i in QI]),
                     'fe': np.array([r['per_q'][i]['f_err'] for i in QI])})
    series[h] = rows


# --- the numbers the assessment rests on ---
print(f'\n{"run":24s} {"blocks":>6s} {"chi2":>9s} {"f0 first->last":>15s} '
      f'{"t(f0=0.5)":>9s} {"ordered":>8s} {"end-to-end":>12s}')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    f0 = np.array([r['f'][0] for r in rows])
    tw = np.array([r['tw'] for r in rows])
    below = np.nonzero(f0 < 0.5)[0]
    if len(below) and below[0] > 0:
        j = below[0]
        t_half = f'{np.interp(0.5, [f0[j], f0[j - 1]], [tw[j], tw[j - 1]]):9.0f}'
    else:
        t_half = '    never'
    arr = [r for r in rows if r['f'][0] < 0.85]            # arrest under way
    ordered = sum(all(r['f'][i] < r['f'][i + 1] for i in range(4)) for r in arr)
    sig = [(r['f'][-1] - r['f'][0]) / np.hypot(r['fe'][-1], r['fe'][0]) for r in arr]
    sig_txt = f'{min(sig):.0f}-{max(sig):.0f} sig' if sig else '-'
    chi = [r['fit']['red_chi2'] for r in rows]
    print(f'{run} {target:5.2f} C  {len(rows):6d} {min(chi):4.2f}-{max(chi):4.2f} '
          f'{f0[0]:6.2f}->{f0[-1]:5.2f}   {t_half} {ordered:4d}/{len(arr):<3d} {sig_txt:>12s}')

print('\nD_med (median time, same per-block fits) before and during conversion:')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    dm = np.array([dmed(r['fit']) for r in rows]) / 1e7
    f0 = np.array([r['f'][0] for r in rows])
    pre, win = dm[f0 > 0.70], dm[f0 <= 0.70]
    txt = f'  {run}: first block {dm[0]:.1f}e7'
    if len(pre) > 1:
        txt += f', last before conversion {pre[-1]:.1f}e7 ({dm[0] / pre[-1]:.1f}x slower)'
    if len(win) > 1:
        txt += (f'; during conversion {win.min():.1f}-{win.max():.1f}e7 ({win.max() / win.min():.1f}x) '
                f'while f0 falls {f0[f0 <= 0.70].max():.2f}->{f0[-1]:.2f}')
    print(txt)

print('\nabsolute, background-subtracted SAXS (mm^-1):')
for fp in sorted(glob.glob(os.path.join(DATA, 'Average_H04*_060C10p_*_results.hdf'))):
    q, I = abs_iq(fp)
    i4, i20 = (np.interp(x, q, I) for x in (0.004, 0.02))
    flat = (q > 0.006) & (q < 0.03)
    print(f'  {os.path.basename(fp)[8:24]} 6 C: I(0.004) {i4:.4f}, I(0.02) {i20:.4f}, '
          f'mean over 0.006-0.03 {np.mean(I[flat]):.4f} +/- {np.std(I[flat]):.4f}')
for h, (run, sample, target) in RUNS.items():
    t, T, I = traces[h]
    hold = t >= 0
    print(f'  {run}: band mean at 20-22 C {np.median(I[T < 22]):.4f}, at the end of the hold '
          f'{np.median(I[hold][-50:]):.2f} ({np.median(I[hold][-50:]) / np.median(I[T < 22]):.0f}x)')


# --- FIGURE ---
apply_style()
fig = plt.figure(figsize=FIG_SIZE)
# Two rows with their own spacing: the top four panels carry their own y labels
# except (c), which shares (b)'s axis; the bottom four share one.
top = GridSpec(1, 4, figure=fig, left=0.075, right=0.985, top=0.955, bottom=0.615,
               wspace=0.62, width_ratios=[1.25, 0.8, 0.8, 1.0])
bot = GridSpec(1, 4, figure=fig, left=0.075, right=0.925, top=0.445, bottom=0.085,
               wspace=0.07)
axa, axb, axc, axd = (fig.add_subplot(top[0, k]) for k in range(4))
axc.sharey(axb)
axf = [fig.add_subplot(bot[0, k]) for k in range(4)]
label_panels([axa, axb, axc, axd] + axf)

# (a) I(Q) through one hold, colour = t_w as in Figure 3a, and the 6 C reference
show = series[SHOW_RUN]
tw_colors = plt.cm.plasma(np.linspace(0.10, 0.88, len(show)))
q, I = abs_iq(REF_6C)
m = I > 0
axa.plot(q[m], I[m], 's', color=COLOR_6C, ms=MS, mfc='none', mew=MEW)
axa.text(0.55, 0.13, '6 °C', transform=axa.transAxes, ha='center', va='bottom',
         color=COLOR_6C)                   # labelled in place, so the key is shorter
handles = []
for r, c in zip(show, tw_colors):
    q, I = abs_iq(r['fp'])
    m = I > 0
    axa.plot(q[m], I[m], 'o', color=c, ms=MS - 1, mfc='none', mew=MEW)
    handles.append(Line2D([], [], marker='o', ls='none', color=c, ms=MS, mfc='none',
                          mew=MEW, label=f'{r["tw"]:.0f} s'))
axa.set_xscale('log')
axa.set_yscale('log')
axa.set_xlim(2.8e-3, 3.6e-2)
q_log_ticks(axa)
axa.set_xlabel(r'$Q$ ($\AA^{-1}$)')
axa.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
add_minor_grid(axa)
lo, hi = axa.get_ylim()
axa.set_ylim(lo, hi * 30)                  # room above the curves for the key
leg = axa.legend(handles=handles, loc='upper right', ncol=3, borderaxespad=0.3,
                 handlelength=0.6, handletextpad=0.25, columnspacing=0.45,
                 labelspacing=0.15)
leg.get_frame().set_linewidth(LW_THIN)
# the temperature sits in the empty band between the 6 C plateau and the hold
axa.text(0.05, 0.29, f'{RUNS[SHOW_RUN][2]:.2f} °C hold', transform=axa.transAxes,
         ha='left', va='center')

# (b) and (c) band intensity on the ramp and through the hold
t_handles = []
for (h, (run, sample, target)), mk in zip(RUNS.items(), T_MARKERS):
    t, T, I = traces[h]
    c = T_COLOR[h]
    ramp = t < 0
    x, y = binned(T[ramp], I[ramp], np.arange(20.0, target + 0.25, 0.25))
    axb.plot(x, y, marker=mk, color=c, ms=MS, mfc='none', mew=MEW, lw=LW_THIN)
    x, y = binned(t[~ramp], I[~ramp], np.arange(0, t.max() + 300, 300))
    axc.plot(x, y, marker=mk, color=c, ms=MS, mfc='none', mew=MEW, lw=LW_THIN)
    t_handles.append(Line2D([], [], marker=mk, color=c, ms=MS, mfc='none', mew=MEW,
                            lw=LW_THIN, label=f'{target:.2f} °C ({sample})'))
axb.set_yscale('log')
axb.set_ylabel(r'$I$, 0.004-0.008 $\AA^{-1}$ (mm$^{-1}$)')
axb.set_xlabel('Zone-3 temperature (°C)')
axb.set_xlim(20.0, 29.0)
axc.set_xlabel(r'$t_w$ (s)')
plt.setp(axc.get_yticklabels(), visible=False)
for ax in (axb, axc):
    add_minor_grid(ax)
leg = axb.legend(handles=t_handles, loc='upper left', borderaxespad=0.3,
                 handlelength=1.4, handletextpad=0.3, labelspacing=0.2)
leg.get_frame().set_linewidth(LW_THIN)

# (d) g2 at the lowest Q through the same hold, colours as in (a)
for r, c in zip(show, tw_colors):
    m = r['tau'] > 0
    axd.errorbar(r['tau'][m], r['g2'][m, 0], yerr=r['g2e'][m, 0], fmt='o', color=c,
                 mfc='none', ms=MS, mew=LW_THIN, capsize=1.0, elinewidth=LW_THIN,
                 alpha=0.75, zorder=2)
    fr, pq = r['fit'], r['fit']['per_q'][0]
    tt = np.logspace(np.log10(r['tau'][m].min()), np.log10(r['tau'][m].max()), 200)
    axd.plot(tt, xpcs_fit.double_exp(tt, pq['tau_fast'], pq['f'], pq['tau_slow'],
                                     fr['p1'], fr['p2']), color=c, lw=LW_DATA, zorder=3)
axd.set_xscale('log')
axd.set_xlim(9e-6, 5.0)
axd.set_ylim(0.99, 1.20)
axd.set_xlabel(r'Delay Time, $\Delta t$ (s)')
axd.set_ylabel(r'$g_2$')
axd.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
axd.xaxis.set_major_formatter(FuncFormatter(
    lambda x, _p: rf'$10^{{{int(round(np.log10(x)))}}}$' if int(round(np.log10(x))) % 2 == 0 else ''))
axd.text(0.96, 0.95, f'$Q$ = {show[0]["q"][0]:.5f} ' + r'$\AA^{-1}$', transform=axd.transAxes,
         ha='right', va='top')
add_minor_grid(axd)

# (e)-(h) f against t_w, marker = Q bin, one panel per temperature
t_max = max(r['tw'] for rows in series.values() for r in rows)
for ax, (h, (run, sample, target)) in zip(axf, RUNS.items()):
    rows = series[h]
    for i, mk in zip(QI, Q_MARKERS):
        x = [r['tw'] for r in rows]
        y = [r['f'][i] for r in rows]
        e = [r['fe'][i] for r in rows]
        ax.plot(x, y, '-', color='0.75', lw=LW_THIN, zorder=1)
        ax.errorbar(x, y, yerr=e, marker=mk, ls='none', color=T_COLOR[h],
                    ms=MS_SPARSE - 1, mfc='none', mew=MEW, capsize=1.0,
                    elinewidth=LW_THIN, zorder=2)
    axr = ax.twinx()
    axr.plot([r['tw'] for r in rows], [dmed(r['fit']) / 1e7 for r in rows], marker='*',
             ls='none', color='k', ms=MS_SPARSE + 1.5, mfc='k', mew=MEW, zorder=4)
    axr.set_ylim(0, D_TOP)
    axr.set_yticks([0, 5, 10, 15])
    if ax is axf[-1]:
        axr.set_ylabel(r'$D_{\mathrm{med}}$ ($10^{7}\ \AA^{2}\,$s$^{-1}$)')
    else:
        axr.set_yticklabels([])
    ax.set_xlim(-400, t_max + 400)
    ax.set_ylim(-0.03, 1.08)
    ax.set_xlabel(r'$t_w$ (s)')
    ax.text(0.04, 0.05, f'{target:.2f} °C', transform=ax.transAxes, ha='left', va='bottom')
    add_minor_grid(ax)
    if ax is not axf[0]:
        ax.set_yticklabels([])
axf[0].set_ylabel('Fast Fraction, $f$')
q0 = series[SHOW_RUN][0]['q']
qh = [Line2D([], [], marker=mk, ls='none', mfc='none', mec='k', mew=MEW,
             ms=MS_SPARSE - 1, label=f'{1e3 * q0[i]:.2f}') for i, mk in zip(QI, Q_MARKERS)]
leg = axf[-1].legend(handles=qh, loc='upper right', ncol=1, borderaxespad=0.3,
                     handlelength=0.8, handletextpad=0.25, labelspacing=0.15,
                     title=r'$Q$ ($10^{-3}\,\AA^{-1}$)')
leg.get_frame().set_linewidth(LW_THIN)

save_fig(fig, 'Candidate_TemperatureSeries_2021-3.pdf')

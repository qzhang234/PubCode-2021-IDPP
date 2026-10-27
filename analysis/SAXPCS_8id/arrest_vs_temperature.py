"""Figure 4: how fast arrest runs, and the single time scale that collapses it.

Three panels, all from the four 2021-3 holds and nothing else.

  (a) the fast fraction f at the largest length scale against t_w, the same
      waiting time Figure 3 uses, one curve per hold;
  (b) the same points with time rescaled by exp[(T - T_REF)/T_s], which brings
      all four onto one master curve;
  (c) the T_s that collapses each Q bin separately, against Q.

WHY f AND NOT THE RATE.  Across these holds the rate of the fast mode stays
inside the same 2-4 x 10^7 A^2/s band once arrest is under way, while the time
for f to reach half its starting value shortens sevenfold over half a degree.
f is also the better determined of the two late in a hold, where the fast decay
is no longer resolved and only the height of the plateau is measured.

WHAT THE COLLAPSE MEANS.  It is a stronger statement than "warmer is faster".
If one rescaling of time superposes the curves, the sample passes through the
SAME sequence of dynamical states at every temperature and temperature sets
only the pace.  That is the dynamic scaling reported for thermoreversible
colloidal gels (Zhang et al., Phys. Rev. Lett. 119, 178006 (2017)), and panel
(c) is the test that matters: one T_s has to work at every length scale, not
just the one plotted.

T_REF IS A CONVENTION, NOT A MEASUREMENT.  Writing the rescaling as
exp[(T - T_REF)/T_s] factorises into exp(-T_REF/T_s) exp(T/T_s), so T_REF
multiplies every series by one common constant: it slides the master curve
sideways and leaves T_s and the collapse quality untouched.  The colloidal-gel
work could anchor it because rheology gave an independent gel temperature.
Nothing here does, so T_REF is set to 28.00 C for readability and no reference
time is quoted from it.

HOW T_s IS FITTED.  One master curve, f = A / (1 + (x/x_c)^n) with
x = t_w exp[(T - T_REF)/T_s), is fitted to the pooled points of all four holds
at once, weighted by the fit uncertainties on f.  T_s is a free parameter of
that single fit alongside A, x_c and n, so it carries a covariance uncertainty
like any other parameter; the error bars drawn in (c) are instead a bootstrap
over the pooled points, which is the larger of the two.

An earlier version scored a trial T_s by binning the rescaled times and
measuring the scatter of f within each bin.  That is NOT used, because the bin
edges were taken from the data range, which itself moves with T_s: points
crossed boundaries discontinuously, the number of contributing bins changed,
and the resulting cost was a step function whose minimum snapped to a handful
of grid values.  It returned identical T_s for neighbouring Q bins and no
uncertainty at all, both of which were artefacts of the binning.

WHAT THE COLLAPSE IS AND IS NOT.  Rescaling improves the master-curve residual
by 2.4 to 3.0 times over no rescaling at every Q, so the effect is real.  The
residual still sits about eight times the fit uncertainties on f, so the
collapse is approximate rather than exact, and T_s drifts monotonically by 13%
across the five Q bins.  The bins share the same correlation functions and
fits, so that drift is not five independent measurements; one T_s describes
them to about 13%, which is the claim the figure supports.

T_s ALSO DEPENDS ON THE CLOCK.  Measured against t_w it is 0.22 K; against
time since the RTD passed 28.00 C it is 0.27 K, with the same collapse
quality.  The slow ramp gives no single unambiguous zero, and that 20% is a
systematic the quoted uncertainty does not cover.

Reads data_2021-3/ only.
"""

import os
import sys
import csv
import glob
from datetime import datetime

import numpy as np
import h5py
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,      # noqa: E402
                              apply_style, add_minor_grid, label_panels, save_fig)
import xpcs_fit                                                           # noqa: E402
from nexus_read import parse_name, read_start_time, read_g2              # noqa: E402

DATA = os.path.join(HERE, 'data_2021-3')
TRACE = os.path.join(DATA, 'trace_2021-3.csv')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')

QI = [0, 1, 2, 3, 4]                  # the five fitted Q bins, as in Figure 3
Q0 = 0                                # the lowest: the largest length scale
ERR_MAX = 0.016                       # noise guard, as in Figure 3
CHI2_MAX = 10.0
AT_TARGET = 0.01                      # C
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
T_REF = 28.00                         # C; see the docstring, it is a convention
TS_SEED = 0.22                        # K, starting guess for the master-curve fit
N_BOOT = 150                          # bootstrap resamples for the T_s error bars

RUNS = {'H0402_D8_280C10p': 28.00, 'H0407_D6_282C10p': 28.25,
        'H0404_D7_285C10p': 28.50, 'H0413_D9_287C10p': 28.75}
T_CMAP, T_LO, T_HI = plt.cm.viridis, 0.05, 0.85
# COLOUR IS TEMPERATURE, MARKER IS Q, everywhere in this figure and in the SI
# companion that repeats (a) and (b) for the four higher Q bins.  Q_MARKERS is
# the same list temperature_series_2021-3.py uses for its five bins, so a shape
# means the same length scale in every figure of the paper.
Q_MARKERS = ['o', 's', '^', 'D', 'v']


# --- beta for this beamtime -------------------------------------------------
# NOT measured over the five fitted bins.  The 2021-3 static reference is the
# NIST glassy carbon standard, whose surface scratches add parasitic scattering
# at low Q.  That scattering is static, so it dilutes the speckle and depresses
# the APPARENT contrast exactly where the paper fits.  The measurement shows it
# plainly: beta climbs from 0.1226 at 0.0036 A^-1 to a plateau of 0.140-0.143
# by 0.0086 A^-1, a 12% rise.  The 2022-1 reference, a different sample, rises
# by only 2% over the same interval, so the effect belongs to the standard and
# not to the beamline.
#
# Contrast is a property of the optics and the detector, near enough constant
# in Q over this range, so it is taken from the clean plateau and applied to
# every bin.  BETA_BAND stops at 0.0165 A^-1, which excludes both a detector
# module gap at 0.0175-0.0185 (beta dips to 0.126) and the stuck pixel's bin at
# 0.0214 (beta reads high).  A single scalar rather than a per-bin value: the
# bin-to-bin scatter inside the plateau is instrumental noise, and letting it
# into f would inject that noise into the very ordering the paper reports.
BETA_BAND = (0.0085, 0.0165)          # A^-1, the flat part of the standard


def flat_level(tau, g, e):
    """Inverse-variance weighted mean of g2 over the measurable delays, and its
    standard error.  Every delay carries its own g2 uncertainty, so each enters
    with weight 1/e^2 and the noisy short delays count for little."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * g[ok]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


def measure_beta(path, label):
    with h5py.File(path, 'r') as hf:
        tau, g2, g2e, q_std = read_g2(hf)
    lv = [flat_level(tau, g2[:, i], g2e[:, i]) for i in range(g2.shape[1])]
    b = np.array([x[0] - 1.0 for x in lv])
    e = np.array([x[1] for x in lv])
    sel = (q_std >= BETA_BAND[0]) & (q_std <= BETA_BAND[1]) & np.isfinite(b)
    w = 1.0 / e[sel] ** 2
    beta = float(np.sum(w * b[sel]) / np.sum(w))
    print(f'beta ({label}) = {beta:.5f} from {sel.sum()} bins over '
          f'{BETA_BAND[0]}-{BETA_BAND[1]} A^-1 (spread {b[sel].min():.4f}-{b[sel].max():.4f})')
    print('   over the five FITTED bins it would read '
          + ', '.join(f'{b[i]:.4f}' for i in QI) + '  <- depressed by the standard')
    return beta


xpcs_fit.CONTRAST = measure_beta(STANDARD, '2021-3 static standard')


def t_zero(run):
    """Start time of the first acquisition at the target: t_w = 0, as in Figure 3."""
    rows = [r for r in csv.DictReader(open(TRACE)) if r['run'] == run]
    rows.sort(key=lambda r: int(r['frame']))
    T = np.array([float(r['qnw3_C']) for r in rows])
    target = float(rows[0]['qnw3_set_C'])
    i = int(np.argmax(T >= target - AT_TARGET))
    return datetime.strptime(rows[i]['start_time'], '%Y-%m-%d %H:%M:%S')


def fit_block(tau, g2, g2e, q):
    """p_fast free, seeded from the minimum of a p_fast profile, as in Figure 3."""
    prof = [(r['red_chi2'], pf, r) for pf in P_GRID
            for r in [xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=pf)]
            if r is not None]
    chi0, pf0, best = min(prof, key=lambda x: x[0])
    xpcs_fit.P_EXP_P0[0] = pf0
    free = xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q)
    xpcs_fit.P_EXP_P0[0] = 0.5
    return free if (free is not None and free['red_chi2'] <= chi0) else best


# --- fit every group of every hold, keeping f in all five Q bins -------------
series, QVALS = {}, None
for run, T in RUNS.items():
    t0 = t_zero(run)
    t, F, E = [], [], []
    for fp in sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')),
                     key=lambda p: parse_name(p)[1]):
        with h5py.File(fp, 'r') as hf:
            w = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if w < 0 or np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue
        r = fit_block(tau, g2, g2e, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        t.append(w)
        F.append([r['per_q'][i]['f'] for i in QI])
        E.append([r['per_q'][i]['f_err'] for i in QI])
        QVALS = [q[i] for i in QI]
    series[run] = (T, np.array(t), np.array(F), np.array(E))
    print(f'  {T:.2f} C: {len(t)} groups, t {t[0]:.0f}-{t[-1]:.0f} s, '
          f'f0 {F[0][Q0]:.2f} -> {F[-1][Q0]:.2f}')


def master_residual(p, T, t, f, e, ts_fixed=None):
    """One master curve f = A/(1+(x/x_c)^n) with x = t_w exp[(T-T_REF)/T_s]."""
    A, lxc, n = p[0], p[1], p[2]
    ts = ts_fixed if ts_fixed is not None else p[3]
    x = t * np.exp((T - T_REF) / ts)
    return (A / (1.0 + (x / np.exp(lxc)) ** n) - f) / e


def pooled(qi):
    T = np.concatenate([np.full(len(t), Tv) for Tv, t, _, _ in series.values()])
    t = np.concatenate([t for _, t, _, _ in series.values()])
    f = np.concatenate([F[:, qi] for _, _, F, _ in series.values()])
    e = np.clip(np.concatenate([E[:, qi] for _, _, _, E in series.values()]), 0.01, None)
    return T, t, f, e


def fit_master(qi, ts_fixed=None, sel=None):
    T, t, f, e = pooled(qi)
    if sel is not None:
        T, t, f, e = T[sel], t[sel], f[sel], e[sel]
    free = ts_fixed is None
    p0 = [1.0, np.log(np.median(t) * 3), 2.0] + ([TS_SEED] if free else [])
    lo = [0.5, 0.0, 0.3] + ([0.05] if free else [])
    hi = [1.3, 20.0, 12.0] + ([3.0] if free else [])
    out = least_squares(master_residual, p0, bounds=(lo, hi), max_nfev=20000,
                        args=(T, t, f, e, ts_fixed))
    rms = float(np.sqrt(np.sum(out.fun ** 2) / max(len(out.fun) - len(p0), 1)))
    return out.x, rms


TS_PER_Q, TS_ERR, RMS_ON, RMS_OFF = [], [], [], []
_rng = np.random.default_rng(0)
for qi in range(len(QI)):
    x, rms = fit_master(qi)
    TS_PER_Q.append(x[3]); RMS_ON.append(rms)
    RMS_OFF.append(fit_master(qi, ts_fixed=1e4)[1])
    n = len(pooled(qi)[0])
    bs = []
    for _ in range(N_BOOT):
        try:
            bs.append(fit_master(qi, sel=_rng.integers(0, n, n))[0][3])
        except Exception:
            pass
    TS_ERR.append(float(np.std(bs, ddof=1)) if len(bs) > 10 else np.nan)
TS_PER_Q = np.array(TS_PER_Q); TS_ERR = np.array(TS_ERR)
TS = float(TS_PER_Q.mean())
TS_SD = float(TS_PER_Q.std(ddof=1))
print('\nT_s from one master curve per Q bin, bootstrap uncertainties:')
for q, s_, e_, a, b in zip(QVALS, TS_PER_Q, TS_ERR, RMS_ON, RMS_OFF):
    print(f'   Q = {q:.5f} A^-1 ({2 * np.pi / q / 10:5.0f} nm): T_s = {s_:.3f} +/- {e_:.3f} K'
          f'   residual {a:.1f} vs {b:.1f} unrescaled ({b / a:.1f}x better)')
print(f'   mean over Q: {TS:.3f} +/- {TS_SD:.3f} K (spread), drift low->high Q '
      f'{100 * (TS_PER_Q[-1] / TS_PER_Q[0] - 1):.0f}%')

# the curve drawn in (b): refitted with T_s held at the across-Q mean
CURVE, _ = fit_master(Q0, ts_fixed=TS)


# --- the 28.00 C hold against the trend set by the other three ---------------
def half_time(t, f):
    thr = f[0] / 2.0
    for i in range(1, len(f)):
        if f[i] <= thr <= f[i - 1]:
            return t[i - 1] + (f[i - 1] - thr) / (f[i - 1] - f[i]) * (t[i] - t[i - 1])
    return np.nan


_arr = [(T, half_time(t, F[:, Q0])) for T, t, F, _ in series.values() if T > T_REF]
_m, _b = np.polyfit([a[0] for a in _arr], [np.log(a[1]) for a in _arr], 1)
_pred = float(np.exp(np.polyval([_m, _b], T_REF)))
_t8 = [(t, F[:, Q0]) for T, t, F, _ in series.values() if T == T_REF][0]
print(f'\nt(f = f0/2) at 28.25/28.50/28.75 C: '
      + ', '.join(f'{a[1]:.0f}' for a in _arr) + f' s  -> T_s = {-1 / _m:.2f} K')
print(f'extrapolated to {T_REF:.2f} C: {_pred:.0f} s; that hold ran to '
      f'{_t8[0][-1]:.0f} s ({100 * _t8[0][-1] / _pred:.0f}% of it) with f still '
      f'{_t8[1][-1]:.2f}, so it is the same curve seen too early')

# --- figures ------------------------------------------------------------------
apply_style()
COLORS = {run: T_CMAP(c) for run, c in zip(RUNS, np.linspace(T_LO, T_HI, len(RUNS)))}


def draw_pair(ax_raw, ax_col, qi, ts):
    """Panels (a) and (b) for one Q bin: measured, then collapsed."""
    mk = Q_MARKERS[qi]
    for run, (T, t, F, E) in series.items():
        c = COLORS[run]
        ax_raw.errorbar(t, F[:, qi], yerr=E[:, qi], fmt=mk, color=c, ms=MS,
                        mfc='none', mew=MEW, lw=LW_THIN, ls='-', elinewidth=LW_THIN,
                        capsize=1.0, capthick=LW_THIN, zorder=3)
        ax_col.errorbar(t * np.exp((T - T_REF) / ts), F[:, qi], yerr=E[:, qi], fmt=mk,
                        color=c, ms=MS, mfc='none', mew=MEW, ls='none',
                        elinewidth=LW_THIN, capsize=1.0, capthick=LW_THIN, zorder=3)
    curve, _ = fit_master(qi, ts_fixed=ts)
    xs = np.logspace(np.log10(min(t.min() * np.exp((T - T_REF) / ts)
                                  for T, t, _, _ in series.values()) * 0.8),
                     np.log10(max(t.max() * np.exp((T - T_REF) / ts)
                                  for T, t, _, _ in series.values()) * 1.2), 300)
    ax_col.plot(xs, curve[0] / (1.0 + (xs / np.exp(curve[1])) ** curve[2]),
                color='0.35', lw=LW_DATA, zorder=2)
    for ax in (ax_raw, ax_col):
        ax.set_ylim(-0.08, 1.15)
        ax.set_ylabel(r'Fast Fraction, $f$')
        add_minor_grid(ax)
    ax_raw.set_xlim(-400, 8900)
    ax_raw.set_xticks([0, 2000, 4000, 6000, 8000])
    ax_raw.set_xlabel(r'$t_w$ (s)')
    ax_col.set_xscale('log')
    ax_col.set_xlabel(r'$t_w\,\mathrm{e}^{(T-T_0)/T_s}$ (s)')


def temperature_key(fig, x_centre, marker):
    """One row: the label, then the four temperatures, colour-coded."""
    h = [Line2D([], [], ls='none', marker='', label=r'$T$ (°C)')]
    h += [Line2D([], [], color=COLORS[run], marker=marker, ls='none', ms=MS,
                 mfc='none', mew=MEW, label=f'{T:.2f}')
          for run, (T, _, _, _) in series.items()]
    leg = fig.legend(handles=h, loc='lower center', bbox_to_anchor=(x_centre, 0.005),
                     ncol=5, handlelength=1.0, columnspacing=1.0, handletextpad=0.4,
                     borderpad=0.3)
    leg.get_frame().set_linewidth(LW_THIN)
    return leg


# ---- Figure 4: the lowest Q bin, plus T_s against Q -------------------------
LEG_H, FIG_H = 0.28, 2.95
fig = plt.figure(figsize=(DOUBLE_COL, FIG_H))
gs = GridSpec(1, 2, figure=fig, left=0.078, right=0.985, top=0.905,
              bottom=0.205 + LEG_H / FIG_H, wspace=0.30)
ax1, ax2 = (fig.add_subplot(gs[0, k]) for k in range(2))
draw_pair(ax1, ax2, Q0, TS)
ax2.text(0.05, 0.06, f'$T_s$ = {TS:.2f} $\\pm$ {TS_SD:.2f} K', transform=ax2.transAxes,
         ha='left', va='bottom')
label_panels([ax1, ax2])
temperature_key(fig, 0.5, Q_MARKERS[Q0])
save_fig(fig, 'Figure4_Arrest_vs_Temperature.pdf')

# ---- SI companion: the same pair for the four higher Q bins -----------------
LEG_H2, FIG_H2 = 0.28, 6.9
fig2 = plt.figure(figsize=(DOUBLE_COL, FIG_H2))
gs2 = GridSpec(3, 4, figure=fig2, left=0.072, right=0.988, top=0.965,
               bottom=0.065 + LEG_H2 / FIG_H2, wspace=0.42, hspace=0.42,
               height_ratios=[1.0, 1.0, 0.78])
axq = []
for k, qi in enumerate([1, 2, 3, 4]):
    a_raw = fig2.add_subplot(gs2[0, k])
    a_col = fig2.add_subplot(gs2[1, k])
    draw_pair(a_raw, a_col, qi, TS_PER_Q[qi])
    a_raw.set_title(f'$Q$ = {QVALS[qi]:.5f} ' + r'$\AA^{-1}$', fontsize=plt.rcParams['font.size'])
    a_col.text(0.05, 0.06, f'$T_s$ = {TS_PER_Q[qi]:.2f} $\\pm$ {TS_ERR[qi]:.2f} K',
               transform=a_col.transAxes, ha='left', va='bottom')
    axq += [a_raw, a_col]
# (i) the scaling temperature that collapses each Q bin on its own
axts = fig2.add_subplot(gs2[2, 1:3])
axts.fill_between([3.2, 8.0], TS - TS_SD, TS + TS_SD, color='crimson', alpha=0.15,
                  lw=0, zorder=0)
axts.axhline(TS, color='crimson', lw=LW_DATA, zorder=1)
for qi, (q, sv, ev) in enumerate(zip(QVALS, TS_PER_Q, TS_ERR)):
    axts.errorbar(q * 1e3, sv, yerr=ev, fmt=Q_MARKERS[qi], color='0.15', ms=MS,
                  mfc='none', mew=MEW, elinewidth=LW_THIN, capsize=1.5,
                  capthick=LW_THIN, ls='none', zorder=3)
axts.set_xlim(3.2, 8.0)
axts.set_ylim(0.10, 0.34)
axts.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
axts.set_ylabel(r'$T_s$ (K)')
add_minor_grid(axts)
label_panels([axq[0], axq[2], axq[4], axq[6], axq[1], axq[3], axq[5], axq[7], axts])
temperature_key(fig2, 0.5, Q_MARKERS[0])
save_fig(fig2, 'FigureS12_Collapse_AllQ.pdf')

"""Figure 3 of the manuscript (file Figure4_Arrest_vs_Temperature.pdf) and its
SI companion, Figure S9 (file FigureS11_Collapse_AllQ.pdf): how fast arrest
runs, the one rescaling of the waiting time that collapses every temperature,
and the SAXS of the same groups on the same clock.

Figure 3, from the four 2021-3 holds and nothing else:
  (a) the fast fraction f at the largest length scale against t_w;
  (b) the same f against X = t_w exp[(T - T0)/T_s], with one master curve
      fitted to all four holds;
  (c) the low-Q SAXS intensity of the same groups against the same X, with T_s
      taken from (b) and not refitted.
SI companion: the four higher Q bins, raw and rescaled with the T_s of (b), the
T_s that collapses each bin on its own, and each bin's time scale X_c at the
T_s of (b).

THE CLOCK.  t_w is the midpoint of the acquisitions a group averages; see
series_2021_3.py.

THE COLLAPSE.  One master curve, f = A / [1 + (X/X_c)^n], with four parameters
(T_s, X_c, n, A), fitted to f at the lowest Q of every group of every hold and
weighted by the fit uncertainties on f floored at 0.01.  T_s carries a
bootstrap uncertainty from resampling whole groups.  Each higher Q bin is
collapsed on its own as well, to show that T_s does not depend on the length
scale; X_c, the time to reach A/2, grows towards higher Q.

T0 IS A CONVENTION.  exp[(T - T0)/T_s] factorises into exp(-T0/T_s) exp(T/T_s),
so T0 multiplies every point by one constant: it slides the master curve along
X and leaves T_s and the collapse untouched.  T0 = 28.00 C, for readability.

THE SAXS ON THE SAME CLOCK.  The mean absolute I(Q) over 0.004-0.008 A^-1 of
each group is placed on t_w exp[(T - T0)/T_s] with the T_s of the f collapse.
The script also finds the T_s that superposes the SAXS on its own (a cubic in
log I against log x, fitted to all four holds), and prints the intensity at
which each hold is half arrested.

WHAT THE COLLAPSE IS NOT.  Rescaling improves the residual about threefold,
but the master curve still sits several fit uncertainties from the points, the
first groups of the 28.50 and 28.75 C holds start below it, and the decline is
broader the warmer the hold.

Reads data_2021-3/ only.
"""

import os
import sys
import glob

import numpy as np
import h5py
from scipy.optimize import least_squares
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,  # noqa: E402
                              apply_style, add_minor_grid, hsv_colors, key_rows,
                              label_panels, save_fig)
import xpcs_fit                                                           # noqa: E402
from nexus_read import parse_name, read_start_time, read_g2              # noqa: E402
from series_2021_3 import (DATA, BAND, T_HUES, T_MARKERS, acquisitions,   # noqa: E402
                           group_midpoint, band_intensity, abs_iq)

STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')

QI = [0, 1, 2, 3, 4]                  # the five fitted Q bins, as in Figure 2
Q0 = 0                                # the lowest: the largest length scale
ERR_MAX = 0.016                       # noise guard, as in Figure 2
CHI2_MAX = 10.0
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
T_REF = 28.00                         # C; a convention, see the docstring
N_BOOT = 200                          # bootstrap resamples

RUNS = {'H0402_D8_280C10p': 28.00, 'H0407_D6_282C10p': 28.25,
        'H0404_D7_285C10p': 28.50, 'H0413_D9_287C10p': 28.75}
# Figure 3 shows one Q, so colour AND marker shape mark temperature there, with
# the hues and shapes the fast-mode figure uses (series_2021_3).  In the SI
# companion each panel is one labelled Q bin, and every point is a circle.

xpcs_fit.CONTRAST = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')


def fit_block(tau, g2, g2e, q):
    """p_fast free, seeded from the minimum of a p_fast profile, as in Figure 2."""
    prof = [(r['red_chi2'], pf, r) for pf in P_GRID
            for r in [xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=pf)]
            if r is not None]
    chi0, pf0, best = min(prof, key=lambda x: x[0])
    xpcs_fit.P_EXP_P0[0] = pf0
    free = xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q)
    xpcs_fit.P_EXP_P0[0] = 0.5
    return free if (free is not None and free['red_chi2'] <= chi0) else best


# --- every group of every hold: f in all five Q bins, and the SAXS ----------
series, saxs, base, QVALS, hold_fps = {}, {}, {}, None, {}
for run, T in RUNS.items():
    frame, t_acq, _, t0 = acquisitions(run)
    files = sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')),
                   key=lambda p: parse_name(p)[1])
    base[run] = band_intensity(files[0])          # the first group of the ramp, ~20 C
    t, F, E, ts_, Is_ = [], [], [], [], []
    hold_fps[run] = []
    for fp in files:
        with h5py.File(fp, 'r') as hf:
            start = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if start < 0:
            continue
        tm = group_midpoint(fp, frame, t_acq)
        ts_.append(tm)
        Is_.append(band_intensity(fp))             # every hold group, fitted or not
        hold_fps[run].append(fp)
        if np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue
        r = fit_block(tau, g2, g2e, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        t.append(tm)
        F.append([r['per_q'][i]['f'] for i in QI])
        E.append([r['per_q'][i]['f_err'] for i in QI])
        QVALS = np.array([q[i] for i in QI])
    series[run] = (T, np.array(t), np.array(F), np.array(E))
    saxs[run] = (T, np.array(ts_), np.array(Is_))
    print(f'  {T:.2f} C: {len(t)} groups, t_w {t[0]:.0f}-{t[-1]:.0f} s, '
          f'f0 {F[0][Q0]:.2f} -> {F[-1][Q0]:.2f}')

# pooled arrays, one row per (group, Q bin)
_T, _t, _q, _f, _e = [], [], [], [], []
for T, t, F, E in series.values():
    for k in range(len(t)):
        for qi in range(len(QI)):
            _T.append(T); _t.append(t[k]); _q.append(QVALS[qi])
            _f.append(F[k, qi]); _e.append(max(E[k, qi], 0.01))
PT, Pt, PQ, PF, PE = map(np.array, (_T, _t, _q, _f, _e))


# --- one master curve per Q bin: its own T_s, or T_s held fixed ---------------
def bin_resid(p, qi, ts_fixed=None, sel=None):
    m = PQ == QVALS[qi]
    T, t, f, e = PT[m], Pt[m], PF[m], PE[m]
    if sel is not None:
        T, t, f, e = T[sel], t[sel], f[sel], e[sel]
    A, lxc, n = p[:3]
    ts = ts_fixed if ts_fixed is not None else p[3]
    return (A / (1.0 + (t * np.exp((T - T_REF) / ts) / np.exp(lxc)) ** n) - f) / e


def fit_bin(qi, ts_fixed=None, sel=None):
    free = ts_fixed is None
    p0 = [0.93, np.log(2e4), 3.5] + ([0.25] if free else [])
    lo = [0.5, 0.0, 0.3] + ([0.05] if free else [])
    hi = [1.3, 25.0, 12.0] + ([3.0] if free else [])
    out = least_squares(bin_resid, p0, bounds=(lo, hi), max_nfev=40000, args=(qi, ts_fixed, sel))
    return out.x, float(np.sum(out.fun ** 2)), len(out.fun)


def bootstrap_ts(qi, rng):
    """T_s of one bin over N_BOOT resamples of its groups (one point per group)."""
    n = int(np.sum(PQ == QVALS[qi]))
    return np.array([fit_bin(qi, sel=rng.integers(0, n, n))[0][3] for _ in range(N_BOOT)])


# Figure 3b: the lowest Q, every hold
_rng = np.random.default_rng(0)
(AMP, LXC, NEXP, TS), CHI_Q0, N_Q0 = fit_bin(Q0)
TS_UNC = float(np.std(bootstrap_ts(Q0, _rng), ddof=1))
CHI_FLAT = fit_bin(Q0, ts_fixed=1e4)[1]          # the same curve without rescaling

# the SI companion: each bin on its own, and each bin's X_c at the T_s of (b)
TS_BIN, TS_BIN_ERR, XC_BIN, XC_BIN_ERR, GAIN = [], [], [], [], []
for qi in range(len(QI)):
    x, chi, n = fit_bin(qi)
    TS_BIN.append(x[3])
    TS_BIN_ERR.append(TS_UNC if qi == Q0 else np.std(bootstrap_ts(qi, _rng), ddof=1))
    GAIN.append(np.sqrt(fit_bin(qi, ts_fixed=1e4)[1] / chi))
    XC_BIN.append(np.exp(fit_bin(qi, ts_fixed=TS)[0][1]))
    bs_xc = [np.exp(fit_bin(qi, ts_fixed=TS, sel=_rng.integers(0, n, n))[0][1]) for _ in range(N_BOOT)]
    XC_BIN_ERR.append(np.std(bs_xc, ddof=1))
TS_BIN, TS_BIN_ERR = np.array(TS_BIN), np.array(TS_BIN_ERR)
XC_BIN, XC_BIN_ERR = np.array(XC_BIN), np.array(XC_BIN_ERR)

print(f'\ncollapse at the lowest Q ({QVALS[Q0]:.5f} A^-1, {2 * np.pi / QVALS[Q0] / 10:.0f} nm), '
      f'all {len(RUNS)} holds ({N_Q0} groups):')
print(f'   T_s = {TS:.3f} +/- {TS_UNC:.3f} K,  X_c = {np.exp(LXC):.0f} s, n = {NEXP:.2f}, A = {AMP:.3f}')
print(f'   residual {np.sqrt(CHI_Q0 / (N_Q0 - 4)):.1f} fit uncertainties, against '
      f'{np.sqrt(CHI_FLAT / (N_Q0 - 3)):.1f} without rescaling '
      f'({np.sqrt(CHI_FLAT / (N_Q0 - 3)) / np.sqrt(CHI_Q0 / (N_Q0 - 4)):.1f}x better)')
print(f'   rate ratio over 0.25 K: {np.exp(0.25 / TS):.1f}x; over 0.5 K: {np.exp(0.5 / TS):.1f}x')
print('\nT_s collapsing each Q bin on its own (bootstrap), and each bin\'s X_c at the T_s of Figure 3b:')
for q, s_, e_, g_, xc, xe in zip(QVALS, TS_BIN, TS_BIN_ERR, GAIN, XC_BIN, XC_BIN_ERR):
    print(f'   Q = {q:.5f} A^-1 ({2 * np.pi / q / 10:4.0f} nm): T_s = {s_:.3f} +/- {e_:.3f} K, '
          f'residual {g_:.1f}x better than unrescaled; X_c = {xc:6.0f} +/- {xe:4.0f} s')
print(f'   per-bin T_s range {TS_BIN.min():.2f}-{TS_BIN.max():.2f} K; '
      f'X_c at {2 * np.pi / QVALS[-1] / 10:.0f} nm is {XC_BIN[-1] / XC_BIN[0]:.2f}x that at '
      f'{2 * np.pi / QVALS[0] / 10:.0f} nm')

# the first group of each hold against the master curve
print('\nfirst group of each hold at the lowest Q, data - master curve:')
for run, (T, t, F, E) in series.items():
    fc = AMP / (1 + (t[0] * np.exp((T - T_REF) / TS) / np.exp(LXC)) ** NEXP)
    print(f'   {T:.2f} C at t_w {t[0]:.0f} s: f {F[0, Q0]:.2f}, curve {fc:.2f} ({F[0, Q0] - fc:+.2f})')
print(f'28.00 C reaches f = 0.5 on this curve at t_w = {np.exp(LXC) * (AMP / 0.5 - 1) ** (1 / NEXP):.0f} s; '
      f'its hold ran to {series["H0402_D8_280C10p"][1][-1]:.0f} s with f still '
      f'{series["H0402_D8_280C10p"][2][-1, Q0]:.2f}')


# --- the SAXS of the same groups, on the same clock --------------------------
def saxs_rss(ts):
    lx = np.concatenate([np.log(t * np.exp((T - T_REF) / ts)) for T, t, _ in saxs.values()])
    ly = np.concatenate([np.log(I) for _, _, I in saxs.values()])
    return float(np.sum((ly - np.polyval(np.polyfit(lx, ly, 3), lx)) ** 2))


_grid = np.linspace(0.10, 1.50, 281)
_rss = np.array([saxs_rss(v) for v in _grid])
TS_SAXS = float(_grid[np.argmin(_rss)])
print(f'\nSAXS band intensity ({BAND[0]}-{BAND[1]} A^-1), all hold groups on t_w exp[(T-T0)/T_s]:')
print(f'   the T_s that superposes it on its own: {TS_SAXS:.2f} K;  scatter about one curve at the '
      f'f value {TS:.2f} K is {saxs_rss(TS) / _rss.min():.2f}x the minimum and '
      f'{saxs_rss(1e4) / saxs_rss(TS):.0f}x smaller than without rescaling')
HALF = {}
for run, (T, t, F, E) in series.items():
    f0 = F[:, Q0]
    j = np.nonzero(f0 < 0.5)[0]
    Ts, ts_s, Is = saxs[run]
    if not len(j):
        print(f'   {T:.2f} C: f0 never reaches 0.5; band {Is[-1]:.1f} mm^-1 at the end (f0 {f0[-1]:.2f}), '
              f'{Is[-1] / base[run]:.0f}x its ~20 C level of {base[run]:.3f}')
        continue
    j = j[0]
    th = np.interp(0.5, [f0[j], f0[j - 1]], [t[j], t[j - 1]])
    Ih = float(np.interp(th, ts_s, Is))
    HALF[T] = (th, Ih)
    print(f'   {T:.2f} C: half arrested at t_w {th:.0f} s, band {Ih:.1f} mm^-1 = {Ih / base[run]:.0f}x its '
          f'~20 C level; then {Is[-1] / Ih:.1f}x more by the end, while f0 falls 0.50 -> {f0[-1]:.2f}')
print(f'   time to f0 = 0.5 shortens {HALF[28.25][0] / HALF[28.75][0]:.1f}-fold from 28.25 to 28.75 C')


# --- the two choices behind (c), and how little rides on them ---------------
# Half arrest is where f at the lowest Q crosses 0.5, interpolated between the
# two groups that bracket it: half of the fluctuation amplitude at 175 nm still
# relaxes on the fast time scale, half is arrested on the 1.6 s window.  Any
# other threshold on the fall of f tells the same story.
print('\nthe half-arrest threshold: band intensity where f0 crosses it, and the growth after it')
for thr in (0.3, 0.4, 0.5, 0.6, 0.7):
    cells = []
    for run, (T, t, F, E) in series.items():
        f0 = F[:, Q0]
        j = np.nonzero(f0 < thr)[0]
        if not len(j) or j[0] == 0:
            continue
        j = j[0]
        th = np.interp(thr, [f0[j], f0[j - 1]], [t[j], t[j - 1]])
        _, ts_s, Is = saxs[run]
        Ih = float(np.interp(th, ts_s, Is))
        cells.append(f'{T:.2f} C {Ih:.1f} mm^-1, then {Is[-1] / Ih:.2f}x')
    print(f'   f0 = {thr:.1f}: ' + '; '.join(cells))

# BAND spans the length scales of the five fitted Q bins, where the intensity
# grows most.  Other bands, inside and beyond them, superpose at the same T_s.
print('\nthe band: T_s that superposes the SAXS on its own, other Q bands')
for lo_, hi_ in [(0.0031, 0.0081), (0.005, 0.010), (0.008, 0.016), (0.012, 0.030)]:
    alt = {}
    for run, (T, ts_s, _) in saxs.items():
        vals = []
        for fp in hold_fps[run]:
            q_, I_ = abs_iq(fp)
            vals.append(float(np.mean(I_[(q_ >= lo_) & (q_ < hi_)])))
        alt[run] = (T, ts_s, np.array(vals))

    def alt_rss(ts):
        lx = np.concatenate([np.log(t * np.exp((T - T_REF) / ts)) for T, t, _ in alt.values()])
        ly = np.concatenate([np.log(I) for _, _, I in alt.values()])
        return float(np.sum((ly - np.polyval(np.polyfit(lx, ly, 3), lx)) ** 2))
    r_ = np.array([alt_rss(v) for v in _grid])
    print(f'   {lo_:.4f}-{hi_:.4f} A^-1: T_s {_grid[np.argmin(r_)]:.2f} K, scatter '
          f'{alt_rss(1e4) / r_.min():.0f}x smaller than without rescaling')


# --- figures ------------------------------------------------------------------
apply_style()
COLORS = dict(zip(RUNS, hsv_colors(T_HUES)))
MARK = dict(zip(RUNS, T_MARKERS))
YLIM_F = (-0.08, 1.15)


def master(xs, A, lxc, n):
    return A / (1.0 + (xs / np.exp(lxc)) ** n)


def raw_panel(ax, qi, marker=None, ms=MS):
    """f against t_w for one Q bin, one line per hold: circles, or the shape
    `marker` gives each run."""
    for run, (T, t, F, E) in series.items():
        ax.errorbar(t, F[:, qi], yerr=E[:, qi], fmt=marker[run] if marker else 'o', color=COLORS[run], ms=ms,
                    mfc='none', mew=MEW, lw=LW_THIN, ls='-', elinewidth=LW_THIN,
                    capsize=1.0, capthick=LW_THIN, zorder=3)
    ax.set_ylim(*YLIM_F)
    ax.set_xlim(-400, 8900)
    ax.set_xticks([0, 2000, 4000, 6000, 8000])
    ax.set_xlabel(r'$t_w$ (s)')
    ax.set_ylabel(r'Fast Fraction, $f$')
    add_minor_grid(ax)


T_TITLE = r'$T$ (°C)'


def temperature_handles(marker=None, ms=MS):
    return [Line2D([], [], color=COLORS[run], marker=marker[run] if marker else 'o', ls='none', ms=ms,
                   mfc='none', mew=MEW, label=f'{T:.2f}') for run, T in RUNS.items()]


# ---- Figure 3 ------------------------------------------------------------------
fig = plt.figure(figsize=(DOUBLE_COL, 3.0))
gs = GridSpec(1, 3, figure=fig, left=0.068, right=0.99, top=0.91, bottom=0.295, wspace=0.36)
ax1, ax2, ax3 = (fig.add_subplot(gs[0, k]) for k in range(3))
raw_panel(ax1, Q0, MARK, MS_SPARSE)

# (b) the largest length scale of every hold on one clock
for run, (T, t, F, E) in series.items():
    ax2.errorbar(t * np.exp((T - T_REF) / TS), F[:, Q0], yerr=E[:, Q0], fmt=MARK[run],
                 color=COLORS[run], ms=MS_SPARSE, mfc='none', mew=MEW, ls='none',
                 elinewidth=LW_THIN, capsize=1.0, capthick=LW_THIN, zorder=3)
xx = np.logspace(3.0, np.log10(8e4), 300)
ax2.plot(xx, master(xx, AMP, LXC, NEXP), color='0.35', lw=LW_DATA, zorder=2)
ax2.set_xscale('log')
ax2.set_xlim(1e3, 8e4)
ax2.set_ylim(*YLIM_F)
ax2.set_xlabel(r'$t_w\,\mathrm{e}^{\Delta T/T_s}$ (s)')
ax2.set_ylabel(r'Fast Fraction, $f$')
ax2.text(0.04, 0.05, f'$T_s$ = {TS:.2f} $\\pm$ {TS_UNC:.2f} K', transform=ax2.transAxes,
         ha='left', va='bottom')
add_minor_grid(ax2)

# (c) the SAXS of the same groups, on the clock of (b)
for run, (T, t, I) in saxs.items():
    ax3.plot(t * np.exp((T - T_REF) / TS), I, MARK[run], color=COLORS[run], ms=MS_SPARSE, mfc='none',
             mew=MEW, zorder=3)
    if T in HALF:
        th, Ih = HALF[T]
        ax3.plot(th * np.exp((T - T_REF) / TS), Ih, '*', color=COLORS[run], mec='k', mew=MEW,
                 ms=MS_SPARSE + 3.5, zorder=4)
ax3.set_xscale('log')
ax3.set_yscale('log')
ax3.set_xlabel(r'$t_w\,\mathrm{e}^{\Delta T/T_s}$ (s)')
ax3.set_ylabel(r'Low-$Q$ Intensity (mm$^{-1}$)')
add_minor_grid(ax3)
label_panels([ax1, ax2, ax3])
star = Line2D([], [], ls='none', marker='*', color='0.7', mec='k', mew=MEW, ms=MS_SPARSE + 3.5,
              label=r'half arrested')
key_rows(fig, [(T_TITLE, temperature_handles(MARK, MS_SPARSE) + [star])], y=0.005)
save_fig(fig, 'Figure4_Arrest_vs_Temperature.pdf')

# ---- SI companion -------------------------------------------------------------
FIG_H2 = 6.9
fig2 = plt.figure(figsize=(DOUBLE_COL, FIG_H2))
gs2 = GridSpec(3, 4, figure=fig2, left=0.072, right=0.988, top=0.965, bottom=0.115,
               wspace=0.42, hspace=0.45, height_ratios=[1.0, 1.0, 0.85])
axes = []
for k, qi in enumerate([1, 2, 3, 4]):
    a_raw = fig2.add_subplot(gs2[0, k])
    a_col = fig2.add_subplot(gs2[1, k])
    raw_panel(a_raw, qi)
    a_raw.set_title(f'$Q$ = {QVALS[qi]:.5f} ' + r'$\AA^{-1}$', fontsize=plt.rcParams['font.size'])
    for run, (T, t, F, E) in series.items():
        a_col.errorbar(t * np.exp((T - T_REF) / TS), F[:, qi], yerr=E[:, qi], fmt='o',
                       color=COLORS[run], ms=MS, mfc='none', mew=MEW, ls='none',
                       elinewidth=LW_THIN, capsize=1.0, capthick=LW_THIN, zorder=3)
    xb, _, _ = fit_bin(qi, ts_fixed=TS)
    xt = np.logspace(2.7, 5.0, 300)
    a_col.plot(xt, master(xt, *xb[:3]), color='0.35', lw=LW_DATA, zorder=2)
    a_col.set_xscale('log')
    a_col.set_xlim(5e2, 1e5)
    a_col.set_ylim(*YLIM_F)
    a_col.set_xlabel(r'$t_w\,\mathrm{e}^{\Delta T/T_s}$ (s)')
    a_col.set_ylabel(r'Fast Fraction, $f$')
    add_minor_grid(a_col)
    axes.append((a_raw, a_col))

# (i) T_s collapsing each bin on its own, against Q; the band is Figure 3b
axts = fig2.add_subplot(gs2[2, 0:2])
axts.fill_between([3.2, 8.0], TS - TS_UNC, TS + TS_UNC, color='crimson', alpha=0.15, lw=0, zorder=0)
axts.axhline(TS, color='crimson', lw=LW_DATA, zorder=1)
for q, sv, ev in zip(QVALS, TS_BIN, TS_BIN_ERR):
    axts.errorbar(q * 1e3, sv, yerr=ev, fmt='o', color='0.15', ms=MS, mfc='none', mew=MEW,
                  elinewidth=LW_THIN, capsize=1.5, capthick=LW_THIN, ls='none', zorder=3)
axts.set_xlim(3.2, 8.0)
axts.set_ylim(0.0, 0.40)
axts.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
axts.set_ylabel(r'$T_s$ (K)')
add_minor_grid(axts)

# (j) each bin's time scale at the T_s of Figure 3b, against Q
axxc = fig2.add_subplot(gs2[2, 2:4])
for q, xc, xe in zip(QVALS, XC_BIN, XC_BIN_ERR):
    axxc.errorbar(q * 1e3, xc / 1e4, yerr=xe / 1e4, fmt='o', color='0.15', ms=MS, mfc='none',
                  mew=MEW, elinewidth=LW_THIN, capsize=1.5, capthick=LW_THIN, ls='none', zorder=3)
axxc.set_xlim(3.2, 8.0)
axxc.set_ylim(0, 3.2)
axxc.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
axxc.set_ylabel(r'$X_c$ ($10^{4}$ s)')
add_minor_grid(axxc)
label_panels([a for a, _ in axes] + [c for _, c in axes] + [axts, axxc])
key_rows(fig2, [(T_TITLE, temperature_handles())], y=0.005)
save_fig(fig2, 'FigureS11_Collapse_AllQ.pdf')

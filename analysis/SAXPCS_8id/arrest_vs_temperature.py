"""Figure 3 (file Figure4_Arrest_vs_Temperature.pdf) and SI Fig. S9 (file
FigureS11_Collapse_AllQ.pdf): how the rate of arrest depends on temperature,
the one rescaling of the waiting time that collapses it, and the low-Q SAXS of
the same groups on the same clock.  Prints the numbers the text quotes.

Figure 3: (a) f at the lowest Q against t_w; (b) the same f against
X = t_w exp[(T - T0)/T_s] with one master curve f = A / [1 + (X/X_c)^n] fitted
to all four isothermal experiments; (c) the low-Q band intensity against X with
the T_s of (b).  SI Fig. S9: the four higher Q bins, raw and rescaled, the T_s
that collapses each bin on its own, and each bin's X_c at the T_s of (b).

T0 = 28.00 C is a convention: it slides the curve along X and leaves T_s
unchanged.  Reads data_2021-3/ only; the README describes the analysis.
"""

import os
import sys

import h5py
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,  # noqa: E402
                              apply_style, add_minor_grid, hsv_colors, key_rows,
                              label_panels, save_fig)
import xpcs_fit                                                           # noqa: E402
from nexus_read import read_start_time, read_g2                           # noqa: E402
from series_2021_3 import (DATA, BAND, T_HUES, T_MARKERS, abs_iq, acquisitions,  # noqa: E402
                           band_intensity, group_files, group_midpoint)

STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')

QI = [0, 1, 2, 3, 4]                  # the five fitted Q bins, as in Figure 2
Q0 = 0                                # the lowest: the largest length scale
ERR_MAX = 0.016                       # noise limit, as in Figure 2
CHI2_MAX = 10.0
T_REF = 28.00                         # C, the convention T0
N_BOOT = 200                          # bootstrap resamples of whole groups
F_ERR_FLOOR = 0.01                    # smallest weight-defining uncertainty on f

RUNS = {'H0402_D8_280C10p': 28.00, 'H0407_D6_282C10p': 28.25,
        'H0404_D7_285C10p': 28.50, 'H0413_D9_287C10p': 28.75}

xpcs_fit.CONTRAST = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')


# --- every group of every isothermal experiment: f in all five bins, and the SAXS
series = {}                           # run -> (T, t_w, f[group, bin], f_err[group, bin]) of the fitted groups
saxs = {}                             # run -> (T, t_w, band intensity) of every group at the target
base = {}                             # run -> band intensity of the first group of the ramp, ~20 C
hold_files = {}                       # run -> every group at the target
D0S = {}                              # run -> D0 (A^2/s) of each fitted group, in the order of series
QVALS = None
for run, T in RUNS.items():
    frame, t_acq, _, t0 = acquisitions(run)
    files = group_files(run)
    base[run] = band_intensity(files[0])
    t_fit = []
    f_fit = []
    e_fit = []
    d_fit = []
    t_all = []
    I_all = []
    hold_files[run] = []
    for fp in files:
        with h5py.File(fp, 'r') as hf:
            start = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if start < 0:
            continue                  # a group of the ramp
        tw = group_midpoint(fp, frame, t_acq)
        t_all.append(tw)
        I_all.append(band_intensity(fp))
        hold_files[run].append(fp)
        if np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue                  # too noisy to fit
        r = xpcs_fit.fit_seeded(tau, g2, g2e, QI, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        t_fit.append(tw)
        f_row = []
        e_row = []
        for i in QI:
            f_row.append(r['per_q'][i]['f'])
            e_row.append(r['per_q'][i]['f_err'])
        f_fit.append(f_row)
        e_fit.append(e_row)
        d_fit.append(r['D0'])
        QVALS = np.array([q[i] for i in QI])
    series[run] = (T, np.array(t_fit), np.array(f_fit), np.array(e_fit))
    D0S[run] = np.array(d_fit)
    saxs[run] = (T, np.array(t_all), np.array(I_all))
    print(f'  {T:.2f} C: {len(t_fit)} groups, t_w {t_fit[0]:.0f}-{t_fit[-1]:.0f} s, '
          f'f0 {f_fit[0][Q0]:.2f} -> {f_fit[-1][Q0]:.2f}')

# every (group, bin) of every isothermal experiment in flat arrays
pool_T = []
pool_t = []
pool_q = []
pool_f = []
pool_e = []
for T, t, F, E in series.values():
    for k in range(len(t)):
        for qi in range(len(QI)):
            pool_T.append(T)
            pool_t.append(t[k])
            pool_q.append(QVALS[qi])
            pool_f.append(F[k, qi])
            pool_e.append(max(E[k, qi], F_ERR_FLOOR))
POOL_T = np.array(pool_T)
POOL_T_W = np.array(pool_t)
POOL_Q = np.array(pool_q)
POOL_F = np.array(pool_f)
POOL_E = np.array(pool_e)


# --- one master curve per Q bin, T_s free or held fixed -----------------------
def bin_resid(p, qi, ts_fixed, sel):
    """Weighted residuals of the master curve in bin qi; p = (A, ln X_c, n[, T_s]).
    sel, if given, picks (and may repeat) the groups used, for the bootstrap."""
    m = POOL_Q == QVALS[qi]
    T = POOL_T[m]
    t = POOL_T_W[m]
    f = POOL_F[m]
    e = POOL_E[m]
    if sel is not None:
        T = T[sel]
        t = t[sel]
        f = f[sel]
        e = e[sel]
    A = p[0]
    lxc = p[1]
    n = p[2]
    if ts_fixed is None:
        ts = p[3]
    else:
        ts = ts_fixed
    return (A / (1.0 + (t * np.exp((T - T_REF) / ts) / np.exp(lxc)) ** n) - f) / e


def fit_bin(qi, ts_fixed=None, sel=None):
    """Master-curve fit of bin qi: (parameters, chi^2, number of points)."""
    if ts_fixed is None:
        p0 = [0.93, np.log(2e4), 3.5, 0.25]
        lo = [0.5, 0.0, 0.3, 0.05]
        hi = [1.3, 25.0, 12.0, 3.0]
    else:
        p0 = [0.93, np.log(2e4), 3.5]
        lo = [0.5, 0.0, 0.3]
        hi = [1.3, 25.0, 12.0]
    out = least_squares(bin_resid, p0, bounds=(lo, hi), max_nfev=40000, args=(qi, ts_fixed, sel))
    return out.x, float(np.sum(out.fun ** 2)), len(out.fun)


def bootstrap_ts(qi, rng):
    """T_s of bin qi over N_BOOT resamples of its groups."""
    n = int(np.sum(POOL_Q == QVALS[qi]))
    values = np.zeros(N_BOOT)
    for b in range(N_BOOT):
        sel = rng.integers(0, n, n)
        values[b] = fit_bin(qi, sel=sel)[0][3]
    return values


# Figure 3b: the lowest Q, every isothermal experiment
rng = np.random.default_rng(0)
params, CHI_Q0, N_Q0 = fit_bin(Q0)
AMP, LXC, NEXP, TS = params
TS_UNC = float(np.std(bootstrap_ts(Q0, rng), ddof=1))
CHI_FLAT = fit_bin(Q0, ts_fixed=1e4)[1]          # the same curve without rescaling

# SI Fig. S9: each bin on its own, and each bin's X_c at the T_s of Figure 3b
TS_BIN = []
TS_BIN_ERR = []
XC_BIN = []
XC_BIN_ERR = []
GAIN = []
for qi in range(len(QI)):
    x, chi, n = fit_bin(qi)
    TS_BIN.append(x[3])
    if qi == Q0:
        TS_BIN_ERR.append(TS_UNC)
    else:
        TS_BIN_ERR.append(np.std(bootstrap_ts(qi, rng), ddof=1))
    GAIN.append(np.sqrt(fit_bin(qi, ts_fixed=1e4)[1] / chi))
    XC_BIN.append(np.exp(fit_bin(qi, ts_fixed=TS)[0][1]))
    xc_resampled = []
    for b in range(N_BOOT):
        sel = rng.integers(0, n, n)
        xc_resampled.append(np.exp(fit_bin(qi, ts_fixed=TS, sel=sel)[0][1]))
    XC_BIN_ERR.append(np.std(xc_resampled, ddof=1))
TS_BIN = np.array(TS_BIN)
TS_BIN_ERR = np.array(TS_BIN_ERR)
XC_BIN = np.array(XC_BIN)
XC_BIN_ERR = np.array(XC_BIN_ERR)

print(f'\ncollapse at the lowest Q ({QVALS[Q0]:.5f} A^-1, {2 * np.pi / QVALS[Q0] / 10:.0f} nm), '
      f'all {len(RUNS)} holds ({N_Q0} groups):')
print(f'   T_s = {TS:.3f} +/- {TS_UNC:.3f} K,  X_c = {np.exp(LXC):.0f} s, n = {NEXP:.2f}, A = {AMP:.3f}')
residual_scaled = np.sqrt(CHI_Q0 / (N_Q0 - 4))
residual_flat = np.sqrt(CHI_FLAT / (N_Q0 - 3))
print(f'   residual {residual_scaled:.1f} fit uncertainties, against '
      f'{residual_flat:.1f} without rescaling '
      f'({residual_flat / residual_scaled:.1f}x better)')
print(f'   rate ratio over 0.25 K: {np.exp(0.25 / TS):.1f}x; over 0.5 K: {np.exp(0.5 / TS):.1f}x')
print('\nT_s collapsing each Q bin on its own (bootstrap), and each bin\'s X_c at the T_s of Figure 3b:')
for k in range(len(QI)):
    q = QVALS[k]
    print(f'   Q = {q:.5f} A^-1 ({2 * np.pi / q / 10:4.0f} nm): T_s = {TS_BIN[k]:.3f} +/- {TS_BIN_ERR[k]:.3f} K, '
          f'residual {GAIN[k]:.1f}x better than unrescaled; X_c = {XC_BIN[k]:6.0f} +/- {XC_BIN_ERR[k]:4.0f} s')
print(f'   per-bin T_s range {TS_BIN.min():.2f}-{TS_BIN.max():.2f} K; '
      f'X_c at {2 * np.pi / QVALS[-1] / 10:.0f} nm is {XC_BIN[-1] / XC_BIN[0]:.2f}x that at '
      f'{2 * np.pi / QVALS[0] / 10:.0f} nm')

print('\nfirst group of each hold at the lowest Q, data - master curve:')
for run, (T, t, F, E) in series.items():
    f_curve = AMP / (1 + (t[0] * np.exp((T - T_REF) / TS) / np.exp(LXC)) ** NEXP)
    print(f'   {T:.2f} C at t_w {t[0]:.0f} s: f {F[0, Q0]:.2f}, curve {f_curve:.2f} ({F[0, Q0] - f_curve:+.2f})')
T_2800, t_2800, F_2800, E_2800 = series['H0402_D8_280C10p']
print(f'28.00 C reaches f = 0.5 on this curve at t_w = {np.exp(LXC) * (AMP / 0.5 - 1) ** (1 / NEXP):.0f} s; '
      f'its hold ran to {t_2800[-1]:.0f} s with f still {F_2800[-1, Q0]:.2f}')


# --- the SAXS of the same groups, on the same clock --------------------------
def saxs_rss(ts, curves):
    """Scatter of log I about one cubic in log X, X = t_w exp[(T - T0)/ts],
    over the curves of every isothermal experiment."""
    lx_parts = []
    ly_parts = []
    for T, t, I in curves.values():
        lx_parts.append(np.log(t * np.exp((T - T_REF) / ts)))
        ly_parts.append(np.log(I))
    lx = np.concatenate(lx_parts)
    ly = np.concatenate(ly_parts)
    return float(np.sum((ly - np.polyval(np.polyfit(lx, ly, 3), lx)) ** 2))


TS_GRID = np.linspace(0.10, 1.50, 281)
rss = np.zeros(len(TS_GRID))
for k in range(len(TS_GRID)):
    rss[k] = saxs_rss(TS_GRID[k], saxs)
TS_SAXS = float(TS_GRID[np.argmin(rss)])
print(f'\nSAXS band intensity ({BAND[0]}-{BAND[1]} A^-1), all hold groups on t_w exp[(T-T0)/T_s]:')
print(f'   the T_s that superposes it on its own: {TS_SAXS:.2f} K;  scatter about one curve at the '
      f'f value {TS:.2f} K is {saxs_rss(TS, saxs) / rss.min():.2f}x the minimum and '
      f'{saxs_rss(1e4, saxs) / saxs_rss(TS, saxs):.0f}x smaller than without rescaling')


def crossing(f0, t, threshold):
    """t_w at which f0 first falls below threshold, interpolated between the two
    groups that bracket it, or None if it never does so after the first group."""
    below = np.nonzero(f0 < threshold)[0]
    if len(below) == 0 or below[0] == 0:
        return None
    j = below[0]
    return np.interp(threshold, [f0[j], f0[j - 1]], [t[j], t[j - 1]])


# Half arrest: f at the lowest Q falls to 0.5.
HALF = {}
for run, (T, t, F, E) in series.items():
    f0 = F[:, Q0]
    T_s, t_s, I_s = saxs[run]
    if not np.any(f0 < 0.5):
        print(f'   {T:.2f} C: f0 never reaches 0.5; band {I_s[-1]:.1f} mm^-1 at the end (f0 {f0[-1]:.2f}), '
              f'{I_s[-1] / base[run]:.0f}x its ~20 C level of {base[run]:.3f}')
        continue
    j = np.nonzero(f0 < 0.5)[0][0]
    t_half = np.interp(0.5, [f0[j], f0[j - 1]], [t[j], t[j - 1]])
    I_half = float(np.interp(t_half, t_s, I_s))
    HALF[T] = (t_half, I_half)
    print(f'   {T:.2f} C: half arrested at t_w {t_half:.0f} s, band {I_half:.1f} mm^-1 = {I_half / base[run]:.0f}x its '
          f'~20 C level; then {I_s[-1] / I_half:.1f}x more by the end, while f0 falls 0.50 -> {f0[-1]:.2f}')
print(f'   time to f0 = 0.5 shortens {HALF[28.25][0] / HALF[28.75][0]:.1f}-fold from 28.25 to 28.75 C')

# Any other threshold on the fall of f gives the same picture.
print('\nthe half-arrest threshold: band intensity where f0 crosses it, and the growth after it')
for threshold in (0.3, 0.4, 0.5, 0.6, 0.7):
    cells = []
    for run, (T, t, F, E) in series.items():
        t_cross = crossing(F[:, Q0], t, threshold)
        if t_cross is None:
            continue
        T_s, t_s, I_s = saxs[run]
        I_cross = float(np.interp(t_cross, t_s, I_s))
        cells.append(f'{T:.2f} C {I_cross:.1f} mm^-1, then {I_s[-1] / I_cross:.2f}x')
    print(f'   f0 = {threshold:.1f}: ' + '; '.join(cells))

# BAND spans the length scales of the five fitted bins; other bands superpose at the same T_s.
print('\nthe band: T_s that superposes the SAXS on its own, other Q bands')
for q_lo, q_hi in [(0.0031, 0.0081), (0.005, 0.010), (0.008, 0.016), (0.012, 0.030)]:
    other = {}
    for run, (T, t_s, I_s) in saxs.items():
        values = []
        for fp in hold_files[run]:
            q, I = abs_iq(fp)
            values.append(float(np.mean(I[(q >= q_lo) & (q < q_hi)])))
        other[run] = (T, t_s, np.array(values))
    rss_other = np.zeros(len(TS_GRID))
    for k in range(len(TS_GRID)):
        rss_other[k] = saxs_rss(TS_GRID[k], other)
    print(f'   {q_lo:.4f}-{q_hi:.4f} A^-1: T_s {TS_GRID[np.argmin(rss_other)]:.2f} K, scatter '
          f'{saxs_rss(1e4, other) / rss_other.min():.0f}x smaller than without rescaling')


# --- activation energy, free energy of activation and bond number -------------
# Over 0.75 K the rescaling exp[(T - T0)/T_s] and an Arrhenius factor
# exp(-E_a/RT) have the same slope d ln(rate)/dT, so E_a = R T^2 / T_s.
R_GAS = 8.314462618                    # J mol^-1 K^-1
K_B = 1.380649e-23                     # J K^-1
H_PLANCK = 6.62607015e-34              # J s
T_MEAN = 273.15 + np.mean(list(RUNS.values()))   # K, mean of the four temperatures
EA = R_GAS * T_MEAN**2 / TS            # J mol^-1
EA_UNC = EA * TS_UNC / TS
print(f'\nactivation energy at T = {T_MEAN:.1f} K: E_a = R T^2/T_s = {EA / (R_GAS * T_MEAN):.0f} '
      f'+/- {EA_UNC / (R_GAS * T_MEAN):.0f} k_B T = {EA / 1e6:.2f} +/- {EA_UNC / 1e6:.2f} MJ/mol')

# Eyring: k = (k_B T/h) exp(-dG/RT).  The Arrhenius slope gives the activation
# enthalpy dH = E_a - RT; reading the time to half arrest as 1/k gives dG.
DH_ACT = EA - R_GAS * T_MEAN
print(f'   activation enthalpy dH = E_a - RT = {DH_ACT / 1e6:.2f} MJ/mol')
for T in sorted(HALF):
    t_half = HALF[T][0]
    T_K = 273.15 + T
    DG_ACT = R_GAS * T_K * np.log(K_B * T_K / H_PLANCK * t_half)
    print(f'   {T:.2f} C: t_half {t_half:.0f} s -> dG = {DG_ACT / (R_GAS * T_K):.1f} k_B T = '
          f'{DG_ACT / 1e3:.0f} kJ/mol; T dS = dH - dG = {(DH_ACT - DG_ACT) / 1e6:.2f} MJ/mol, '
          f'dS = {(DH_ACT - DG_ACT) / T_K / 1e3:.1f} kJ/(mol K); entropy offsets '
          f'{100 * (1 - DG_ACT / DH_ACT):.1f}% of dH')
print(f'   a 10^6-fold error in the rate moves dG by {np.log(1e6):.1f} k_B T')

# E_a belongs to the clock (how fast the state evolves), not to the relaxation at
# a given state: compare the fast-mode rate D0 of the four experiments at the
# same rescaled time X.  Only groups with f0 >= 0.25, where the fast decay fixes
# a rate (as in Figure 2c); D0 interpolated in log-log between bracketing groups.
print('\nD0 (1e7 A^2/s) against X = t_w exp[(T - T0)/T_s], groups with f0 >= 0.25:')
X_D0 = {}
for run, (T, t, F, E) in series.items():
    X = t * np.exp((T - T_REF) / TS)
    keep = F[:, Q0] >= 0.25
    X_D0[run] = (T, X[keep], D0S[run][keep])
    cells = []
    for k in np.nonzero(keep)[0]:
        cells.append(f'X {X[k]:.0f} s (f0 {F[k, Q0]:.2f}) {D0S[run][k] / 1e7:.2f}')
    print(f'   {T:.2f} C: ' + '; '.join(cells))
print('   at the same X: D0, the fast relaxation time 1/(D0 Q^2) at the lowest Q, and the waiting')
print('   time t_w = X exp[-(T - T0)/T_s] each temperature needs to reach that stage')
for x_common in (3000.0, 6600.0):
    cells = []
    temps = []
    values = []
    waits = []
    for run, (T, Xk, Dk) in X_D0.items():
        if len(Xk) < 2 or x_common < Xk.min() or x_common > Xk.max():
            continue
        d = np.exp(np.interp(np.log(x_common), np.log(Xk), np.log(Dk)))
        tau = 1.0 / (d * QVALS[Q0]**2)            # s, the 1/e time of the fast decay
        t_wait = x_common * np.exp(-(T - T_REF) / TS)
        cells.append(f'{T:.2f} C: D0 {d / 1e7:.2f}, 1/(D0 Q^2) {tau * 1e3:.2f} ms, t_w {t_wait:.0f} s')
        temps.append(T)
        values.append(d)
        waits.append(t_wait)
    spread = max(values) / min(values)
    clock = max(waits) / min(waits)
    print(f'   X = {x_common:.0f} s: ' + '; '.join(cells))
    print(f'      relaxation time spread {spread:.1f}x, waiting time {clock:.1f}x over the same temperatures')

# Bond number N = E_a / (transition enthalpy of one pentapeptide).
# Calorimetry exists only for long poly(VPAVG) chains in water: J/g on heating,
# and the number of repeats per chain.
M_REPEAT = 99.131 + 97.115 + 71.078 + 99.131 + 57.051   # g/mol, residues V P A V G
CAL = 4.184                                            # J per cal
DH_LONG = [('Luan & Urry 1991, synthetic, >50 kDa', 6.52 * CAL, 118),
           ('Reguera 2003, synthetic, Mn 74.8 kDa', 34.0, 74777 / M_REPEAT),
           ('Ribeiro 2009, recombinant (VPAVG)62', 34.0, 62)]
# How the transition enthalpy of a related ELP, [(VPGVG)2 VPGEG (VPGVG)2]n at
# pH 2.5, depends on chain length (Girotti 2004, Fig. 4, read off the plot):
# pentapeptides per chain, and J/g.  Linear interpolation between its points
# scales each long-chain value to 30 repeats.
GIROTTI_N = np.array([25, 45, 75, 150, 225])
GIROTTI_DH = np.array([4.49, 7.82, 8.24, 8.87, 9.03])
N_REP = 30
print(f'\nbond number N = E_a / dH per pentapeptide ({M_REPEAT:.2f} g/mol), chains = N / {N_REP}:')
for source, dh_g, n_chain in DH_LONG:
    scale = np.interp(N_REP, GIROTTI_N, GIROTTI_DH) / np.interp(n_chain, GIROTTI_N, GIROTTI_DH)
    dh_long = dh_g * M_REPEAT              # J/mol per pentapeptide, long chain
    dh_30 = dh_long * scale                # J/mol per pentapeptide, scaled to 30 repeats
    print(f'   {source} ({n_chain:.0f} repeats): {dh_g:.1f} J/g = {dh_long / 1e3:.2f} kJ/mol '
          f'({dh_long / (R_GAS * T_MEAN):.2f} k_B T) -> N {EA / dh_long:.0f}, {EA / dh_long / N_REP:.1f} chains; '
          f'scaled to {N_REP} repeats (x{scale:.3f}): {dh_g * scale:.1f} J/g = {dh_30 / 1e3:.2f} kJ/mol -> '
          f'N {EA / dh_30:.0f}, {EA / dh_30 / N_REP:.1f} chains')
print(f'   Luan & Urry 1991: 0.5 M NaCl raises the enthalpy {3.01 / 2.76:.3f}-fold (2.76 -> 3.01 kcal/mol)')
print(f'   the desolvation transition (Reguera 2003, Friedman E_a 220-480 kJ/mol) is '
      f'{EA / 480e3:.1f}-{EA / 220e3:.1f}x less steep; its T_s would be '
      f'{R_GAS * T_MEAN**2 / 480e3:.1f}-{R_GAS * T_MEAN**2 / 220e3:.1f} K')


# --- figures ------------------------------------------------------------------
apply_style()
COLORS = {}
MARK = {}
colors = hsv_colors(T_HUES)
run_names = list(RUNS)
for k in range(len(run_names)):
    COLORS[run_names[k]] = colors[k]
    MARK[run_names[k]] = T_MARKERS[k]
YLIM_F = (-0.08, 1.15)
T_TITLE = r'$T$ (°C)'


def master(xs, A, lxc, n):
    return A / (1.0 + (xs / np.exp(lxc)) ** n)


def raw_panel(ax, qi, marker=None, ms=MS):
    """f against t_w in bin qi, one line per temperature: circles, or the shape
    marker[run] gives each run."""
    for run, (T, t, F, E) in series.items():
        if marker is None:
            fmt = 'o'
        else:
            fmt = marker[run]
        ax.errorbar(t, F[:, qi], yerr=E[:, qi], fmt=fmt, color=COLORS[run], ms=ms,
                    mfc='none', mew=MEW, lw=LW_THIN, ls='-', elinewidth=LW_THIN,
                    capsize=1.0, capthick=LW_THIN, zorder=3)
    ax.set_ylim(*YLIM_F)
    ax.set_xlim(-400, 8900)
    ax.set_xticks([0, 2000, 4000, 6000, 8000])
    ax.set_xlabel(r'$t_w$ (s)')
    ax.set_ylabel(r'Fast Fraction, $f$')
    add_minor_grid(ax)


def temperature_handles(marker=None, ms=MS):
    handles = []
    for run, T in RUNS.items():
        if marker is None:
            mk = 'o'
        else:
            mk = marker[run]
        handles.append(Line2D([], [], color=COLORS[run], marker=mk, ls='none', ms=ms,
                              mfc='none', mew=MEW, label=f'{T:.2f}'))
    return handles


# ---- Figure 3 ------------------------------------------------------------------
fig = plt.figure(figsize=(DOUBLE_COL, 3.0))
gs = GridSpec(1, 3, figure=fig, left=0.068, right=0.99, top=0.91, bottom=0.295, wspace=0.36)
ax1 = fig.add_subplot(gs[0, 0])
ax2 = fig.add_subplot(gs[0, 1])
ax3 = fig.add_subplot(gs[0, 2])
raw_panel(ax1, Q0, MARK, MS_SPARSE)

# (b) the largest length scale of every isothermal experiment on one clock
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
        t_half, I_half = HALF[T]
        ax3.plot(t_half * np.exp((T - T_REF) / TS), I_half, '*', color=COLORS[run], mec='k', mew=MEW,
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
save_fig(fig, os.path.join(HERE, 'Figure4_Arrest_vs_Temperature.pdf'))

# ---- SI Fig. S9 ----------------------------------------------------------------
fig2 = plt.figure(figsize=(DOUBLE_COL, 6.9))
gs2 = GridSpec(3, 4, figure=fig2, left=0.072, right=0.988, top=0.965, bottom=0.115,
               wspace=0.42, hspace=0.45, height_ratios=[1.0, 1.0, 0.85])
raw_axes = []
collapsed_axes = []
for k in range(4):
    qi = k + 1                        # the four higher Q bins
    a_raw = fig2.add_subplot(gs2[0, k])
    a_col = fig2.add_subplot(gs2[1, k])
    raw_panel(a_raw, qi)
    a_raw.set_title(f'$Q$ = {QVALS[qi]:.5f} ' + r'$\AA^{-1}$', fontsize=plt.rcParams['font.size'])
    for run, (T, t, F, E) in series.items():
        a_col.errorbar(t * np.exp((T - T_REF) / TS), F[:, qi], yerr=E[:, qi], fmt='o',
                       color=COLORS[run], ms=MS, mfc='none', mew=MEW, ls='none',
                       elinewidth=LW_THIN, capsize=1.0, capthick=LW_THIN, zorder=3)
    xb = fit_bin(qi, ts_fixed=TS)[0]
    xt = np.logspace(2.7, 5.0, 300)
    a_col.plot(xt, master(xt, xb[0], xb[1], xb[2]), color='0.35', lw=LW_DATA, zorder=2)
    a_col.set_xscale('log')
    a_col.set_xlim(5e2, 1e5)
    a_col.set_ylim(*YLIM_F)
    a_col.set_xlabel(r'$t_w\,\mathrm{e}^{\Delta T/T_s}$ (s)')
    a_col.set_ylabel(r'Fast Fraction, $f$')
    add_minor_grid(a_col)
    raw_axes.append(a_raw)
    collapsed_axes.append(a_col)

# (i) T_s collapsing each bin on its own, against Q; the line and band are Figure 3b
axts = fig2.add_subplot(gs2[2, 0:2])
axts.fill_between([3.2, 8.0], TS - TS_UNC, TS + TS_UNC, color='crimson', alpha=0.15, lw=0, zorder=0)
axts.axhline(TS, color='crimson', lw=LW_DATA, zorder=1)
for k in range(len(QI)):
    axts.errorbar(QVALS[k] * 1e3, TS_BIN[k], yerr=TS_BIN_ERR[k], fmt='o', color='0.15', ms=MS,
                  mfc='none', mew=MEW, elinewidth=LW_THIN, capsize=1.5, capthick=LW_THIN,
                  ls='none', zorder=3)
axts.set_xlim(3.2, 8.0)
axts.set_ylim(0.0, 0.40)
axts.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
axts.set_ylabel(r'$T_s$ (K)')
add_minor_grid(axts)

# (j) each bin's time scale at the T_s of Figure 3b, against Q
axxc = fig2.add_subplot(gs2[2, 2:4])
for k in range(len(QI)):
    axxc.errorbar(QVALS[k] * 1e3, XC_BIN[k] / 1e4, yerr=XC_BIN_ERR[k] / 1e4, fmt='o', color='0.15',
                  ms=MS, mfc='none', mew=MEW, elinewidth=LW_THIN, capsize=1.5, capthick=LW_THIN,
                  ls='none', zorder=3)
axxc.set_xlim(3.2, 8.0)
axxc.set_ylim(0, 3.2)
axxc.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
axxc.set_ylabel(r'$X_c$ ($10^{4}$ s)')
add_minor_grid(axxc)
label_panels(raw_axes + collapsed_axes + [axts, axxc])
key_rows(fig2, [(T_TITLE, temperature_handles())], y=0.005)
save_fig(fig2, os.path.join(HERE, 'FigureS11_Collapse_AllQ.pdf'))

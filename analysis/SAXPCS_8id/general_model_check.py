"""Does the collapse of Figure 3b survive the choices in the g2 fit?

Prints the numbers SI sections 4.3 and 6.2 quote; draws nothing; takes about
six minutes.  The collapse (one master curve through f at the lowest Q of every
group of the four isothermal experiments) is redone with f from five versions
of the fit:

  diff      Equation 2 as every figure uses it: tau_fast = 1/(D0 Q^2), p_fast
            free and seeded from a p_fast profile;
  seeded    the general two-mode form, tau_fast free in every Q bin, seeded
            the same way;
  default   the general form started from the bare default p_fast = 0.5;
  p05       Equation 2 with p_fast fixed at 0.5;
  beta_low  Equation 2 at the contrast the standard reads over the five fitted
            bins instead of its plateau.

Reads data_2021-3/ only.
"""

import os
import sys

import h5py
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import xpcs_fit                                                    # noqa: E402
from nexus_read import read_start_time, read_g2                    # noqa: E402
from series_2021_3 import DATA, acquisitions, group_files, group_midpoint  # noqa: E402

STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI = [0, 1, 2, 3, 4]
Q0 = 0
ERR_MAX = 0.016
CHI2_MAX = 10.0
T_REF = 28.00
RUNS = {'H0402_D8_280C10p': 28.00, 'H0407_D6_282C10p': 28.25,
        'H0404_D7_285C10p': 28.50, 'H0413_D9_287C10p': 28.75}

# The contrast every figure uses (the plateau of the static standard), and the
# lower value the same standard reads over the five fitted bins.
BETA = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')
with h5py.File(STANDARD, 'r') as hf:
    tau_std, g2_std, g2e_std, q_std = read_g2(hf)
levels = np.zeros(len(QI))
weights = np.zeros(len(QI))
for k in range(len(QI)):
    level, level_err = xpcs_fit.flat_level(tau_std, g2_std[:, QI[k]], g2e_std[:, QI[k]])
    levels[k] = level - 1
    weights[k] = 1 / level_err ** 2
BETA_LOW = float(np.sum(weights * levels) / np.sum(weights))
print(f'beta over the fitted bins: {BETA_LOW:.4f}, {100 * (BETA_LOW / BETA - 1):+.1f}% against {BETA:.4f}')


def fits(tau, g2, g2e, q):
    """The five versions of the fit of one group, or None if Equation 2 itself
    fails the chi^2 limit."""
    xpcs_fit.CONTRAST = BETA
    out = {}
    out['diff'] = xpcs_fit.fit_seeded(tau, g2, g2e, QI, q, model='diffusive')
    if out['diff'] is None or out['diff']['red_chi2'] >= CHI2_MAX:
        return None
    out['seeded'] = xpcs_fit.fit_seeded(tau, g2, g2e, QI, q, model='general')
    out['default'] = xpcs_fit.fit_g2_global(tau, g2, g2e, QI)        # p_fast starts at 0.5
    out['p05'] = xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=0.5)
    xpcs_fit.CONTRAST = BETA_LOW
    out['beta_low'] = xpcs_fit.fit_seeded(tau, g2, g2e, QI, q, model='diffusive')
    xpcs_fit.CONTRAST = BETA
    return out


KEYS = ('diff', 'seeded', 'default', 'p05', 'beta_low')
S = {}                                # version -> one dict per group
for key in KEYS:
    S[key] = []
QV = None
bad = []                              # groups where the default start lands in a false minimum
for run, T in RUNS.items():
    frame, t_acq, _, t0 = acquisitions(run)
    for fp in group_files(run):
        with h5py.File(fp, 'r') as hf:
            start = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if start < 0 or np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue
        out = fits(tau, g2, g2e, q)
        if out is None:
            continue
        complete = True
        for key in KEYS:
            if out[key] is None:
                complete = False
        if not complete:
            continue
        tw = group_midpoint(fp, frame, t_acq)
        if out['default']['red_chi2'] > 10:
            bad.append((run, tw, out['default']['red_chi2'], out['seeded']['red_chi2']))
        for key in KEYS:
            per_q = out[key]['per_q']
            S[key].append({'T': T, 'tw': tw,
                           'f': [per_q[i]['f'] for i in QI],
                           'f_err': [per_q[i]['f_err'] for i in QI],
                           'tau_fast': [per_q[i]['tau_fast'] for i in QI]})
        QV = np.array([q[i] for i in QI])


def collapse(rows):
    """T_s of the master curve of Figure 3b through f at the lowest Q of rows."""
    T = np.array([r['T'] for r in rows])
    t = np.array([r['tw'] for r in rows])
    f = np.array([r['f'][Q0] for r in rows])
    e = np.clip([r['f_err'][Q0] for r in rows], 0.01, None)

    def residual(p):              # p = (T_s, ln X_c, n, A)
        x = t * np.exp((T - T_REF) / p[0])
        return (p[3] / (1 + (x / np.exp(p[1])) ** p[2]) - f) / e

    out = least_squares(residual, [0.25, np.log(2e4), 3.5, 0.93], max_nfev=40000,
                        bounds=([0.05, 0.0, 0.3, 0.5], [3.0, 25.0, 12.0, 1.3]))
    return out.x[0]


def lowest_q_f(rows):
    return np.array([r['f'][Q0] for r in rows])


def count_rising(rows):
    """Number of groups whose f increases monotonically with Q."""
    n = 0
    for r in rows:
        if np.all(np.diff(r['f']) > 0):
            n += 1
    return n


print(f'\n@@ collapse at the lowest Q ({len(S["diff"])} groups) with f from each version of the g2 fit:')
ref = collapse(S['diff'])
for key in KEYS:
    ts = collapse(S[key])
    print(f'@@  {key:8s}: T_s = {ts:.3f} K ({100 * (ts / ref - 1):+.1f}%)')

f_ref = lowest_q_f(S['diff'])
print('\n@@ |df| at the lowest Q against Equation 2 as fitted')
for key in ('seeded', 'default', 'p05', 'beta_low'):
    d = np.abs(lowest_q_f(S[key]) - f_ref)
    print(f'@@  {key:8s}: at most {d.max():.3f}, median {np.median(d):.3f}'
          f'   (where the fast mode has nearly gone, f < 0.15: at most {d[f_ref < 0.15].max():.3f})')
print(f'@@  f at the lowest Q spans {f_ref.min():.2f}-{f_ref.max():.2f} under Equation 2')

# the arrested fraction at the end of the 28.50 C experiment, and the ordering in Q, at the lower contrast
last = None
for k in range(len(S['diff'])):
    if S['diff'][k]['T'] == 28.50:
        last = k
print(f'@@ beta_low: 1 - f at the end of the 28.50 C hold {1 - S["diff"][last]["f"][Q0]:.2f} -> '
      f'{1 - S["beta_low"][last]["f"][Q0]:.2f}; groups with f rising monotonically in Q: '
      f'{count_rising(S["diff"])} -> {count_rising(S["beta_low"])}'
      f' of {len(S["diff"])}')

# SI section 4.3: the Q exponent of the released tau_fast, from the seeded general fits
log_q = np.log(QV)
slopes = []
slope_f0 = []
for k in range(len(S['seeded'])):
    tau_fast = S['seeded'][k]['tau_fast']
    if np.all(np.isfinite(tau_fast)) and np.all(np.asarray(tau_fast) > 0):
        slopes.append(np.polyfit(log_q, np.log(tau_fast), 1)[0])
        slope_f0.append(S['diff'][k]['f'][Q0])
slopes = np.array(slopes)
slope_f0 = np.array(slope_f0)
resolved = slope_f0 > 0.15
print(f'\n@@ Q exponent of the released tau_fast (seeded): median {np.median(slopes[resolved]):.2f} over the '
      f'{resolved.sum()} groups with f > 0.15 at the lowest Q, 80% between {np.percentile(slopes[resolved], 10):.1f} '
      f'and {np.percentile(slopes[resolved], 90):.1f}')
print(f'\n@@ groups where the DEFAULT seed lands above chi2 10 ({len(bad)}):')
for run, tw, chi_default, chi_seeded in bad:
    print(f'@@   {run} t_w={tw:.0f}s  default chi2 {chi_default:.1f}  vs seeded {chi_seeded:.2f}')

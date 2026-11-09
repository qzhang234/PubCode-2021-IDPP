"""Does the collapse of Figure 3b survive the choices in the g2 fit?

Prints the numbers SI sections 4.3 and 6.2 quote.  The collapse of Figure 3b
-- one master curve through f at the lowest Q of all four holds, with scaling
temperature T_s, on group-midpoint t_w -- is redone with f from five versions
of the g2 fit:

  diff      Equation 2, the fast time tied to 1/(D0 Q^2), p_fast free and
            seeded from a p_fast profile (what every figure uses);
  seeded    the general two-mode form of Bahadur 2019, tau_fast free in every
            Q bin, seeded the same way;
  default   the same general form started from the bare default p_fast = 0.5;
  p05       Equation 2 with p_fast frozen at 0.5 in every group;
  beta_low  Equation 2 at the contrast the standard reads over the five fitted
            bins, where its scratches depress it, instead of its plateau.

WHY THE SEED IS PART OF THE ANSWER.  With p_fast free, the two-mode
likelihood has local minima in which the fast term imitates the slow one.
Every fit in this paper therefore starts from the minimum of a p_fast profile
and is then released, and the released fit is kept only if it does not make
chi2 worse.  Started from the bare default instead, a few groups land in one
of those minima; both versions are printed so the difference cannot be
mistaken for a discrepancy.

Reads data_2021-3/ only.  Draws nothing.  Takes about a quarter of an hour.
"""

import glob
import os
import sys

import h5py
import numpy as np
from scipy.optimize import least_squares

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
import xpcs_fit                                                    # noqa: E402
from nexus_read import parse_name, read_start_time, read_g2        # noqa: E402
from series_2021_3 import DATA, acquisitions, group_midpoint       # noqa: E402

STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI, Q0, ERR_MAX, CHI2_MAX = [0, 1, 2, 3, 4], 0, 0.016, 10.0
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
T_REF = 28.00
RUNS = {'H0402_D8_280C10p': 28.00, 'H0407_D6_282C10p': 28.25,
        'H0404_D7_285C10p': 28.50, 'H0413_D9_287C10p': 28.75}

# The contrast every figure uses, from the plateau of the static standard, and
# the lower value the same standard reads over the five fitted bins.
BETA = xpcs_fit.measure_beta(STANDARD, '2021-3 static standard')
with h5py.File(STANDARD, 'r') as _hf:
    _tau, _g2, _g2e, _ = read_g2(_hf)
_lv = [xpcs_fit.flat_level(_tau, _g2[:, i], _g2e[:, i]) for i in QI]
_w = np.array([1 / e ** 2 for _, e in _lv])
BETA_LOW = float(np.sum(_w * np.array([v - 1 for v, _ in _lv])) / np.sum(_w))
print(f'beta over the fitted bins: {BETA_LOW:.4f}, {100 * (BETA_LOW / BETA - 1):+.1f}% against {BETA:.4f}')


def profiled(fit, *args):
    """Fit started from the minimum of a p_fast profile, then released."""
    prof = [(r['red_chi2'], pf, r) for pf in P_GRID
            for r in [fit(*args, p1_fixed=pf)] if r is not None]
    if not prof:
        return None
    chi0, pf0, best = min(prof, key=lambda x: x[0])
    xpcs_fit.P_EXP_P0[0] = pf0
    free = fit(*args)
    xpcs_fit.P_EXP_P0[0] = 0.5
    return free if (free is not None and free['red_chi2'] <= chi0) else best


def diffusive(tau, g2, g2e, q, p1_fixed=None):
    return xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=p1_fixed)


def general(tau, g2, g2e, q, p1_fixed=None):
    return xpcs_fit.fit_g2_global(tau, g2, g2e, QI, p1_fixed=p1_fixed)


def fits(tau, g2, g2e, q):
    xpcs_fit.CONTRAST = BETA
    out = {'diff': profiled(diffusive, tau, g2, g2e, q)}
    if out['diff'] is None or out['diff']['red_chi2'] >= CHI2_MAX:
        return None
    out['seeded'] = profiled(general, tau, g2, g2e, q)
    xpcs_fit.P_EXP_P0[0] = 0.5
    out['default'] = general(tau, g2, g2e, q)
    out['p05'] = diffusive(tau, g2, g2e, q, p1_fixed=0.5)
    xpcs_fit.CONTRAST = BETA_LOW
    out['beta_low'] = profiled(diffusive, tau, g2, g2e, q)
    xpcs_fit.CONTRAST = BETA
    return out


KEYS = ('diff', 'seeded', 'default', 'p05', 'beta_low')
S = {k: [] for k in KEYS}          # rows of (T, t_w, f[5], f_err[5], tau_fast[5])
QV, bad = None, []
for run, T in RUNS.items():
    frame, t_acq, _, t0 = acquisitions(run)
    for fp in sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')),
                     key=lambda p: parse_name(p)[1]):
        with h5py.File(fp, 'r') as hf:
            start = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        if start < 0 or np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX:
            continue
        out = fits(tau, g2, g2e, q)
        if out is None or any(out[k] is None for k in KEYS):
            continue
        tw = group_midpoint(fp, frame, t_acq)
        if out['default']['red_chi2'] > 10:
            bad.append((run, tw, out['default']['red_chi2'], out['seeded']['red_chi2']))
        for k in KEYS:
            r = out[k]
            S[k].append((T, tw, [r['per_q'][i]['f'] for i in QI], [r['per_q'][i]['f_err'] for i in QI],
                         [r['per_q'][i]['tau_fast'] for i in QI]))
        QV = np.array([q[i] for i in QI])


def collapse(rows):
    """The master curve of Figure 3b, f at the lowest Q of every hold: returns T_s."""
    T = np.array([r[0] for r in rows]); t = np.array([r[1] for r in rows])
    f = np.array([r[2][Q0] for r in rows]); e = np.clip([r[3][Q0] for r in rows], 0.01, None)

    def res(p):
        x = t * np.exp((T - T_REF) / p[0])
        return (p[3] / (1 + (x / np.exp(p[1])) ** p[2]) - f) / e
    out = least_squares(res, [0.25, np.log(2e4), 3.5, 0.93], max_nfev=40000,
                        bounds=([0.05, 0.0, 0.3, 0.5], [3.0, 25.0, 12.0, 1.3]))
    return out.x[0]


print(f'\n@@ collapse at the lowest Q ({len(S["diff"])} groups) with f from each version of the g2 fit:')
ref = collapse(S['diff'])
for k in KEYS:
    ts = collapse(S[k])
    print(f'@@  {k:8s}: T_s = {ts:.3f} K ({100 * (ts / ref - 1):+.1f}%)')

f_ref = np.array([r[2][Q0] for r in S['diff']])
print('\n@@ |df| at the lowest Q against Equation 2 as fitted')
for k in ('seeded', 'default', 'p05', 'beta_low'):
    d = np.abs(np.array([r[2][Q0] for r in S[k]]) - f_ref)
    print(f'@@  {k:8s}: at most {d.max():.3f}, median {np.median(d):.3f}'
          f'   (where the fast mode has nearly gone, f < 0.15: at most {d[f_ref < 0.15].max():.3f})')
print(f'@@  f at the lowest Q spans {f_ref.min():.2f}-{f_ref.max():.2f} under Equation 2')

# the arrested fraction at the end of the 28.50 C hold, and the ordering in Q, at the lower contrast
last = [i for i, r in enumerate(S['diff']) if r[0] == 28.50][-1]
print(f'@@ beta_low: 1 - f at the end of the 28.50 C hold {1 - S["diff"][last][2][Q0]:.2f} -> '
      f'{1 - S["beta_low"][last][2][Q0]:.2f}; groups with f rising monotonically in Q: '
      f'{sum(all(np.diff(r[2]) > 0) for r in S["diff"])} -> {sum(all(np.diff(r[2]) > 0) for r in S["beta_low"])}'
      f' of {len(S["diff"])}')

# SI section 4.3 reports the Q exponent the released tau_fast returns, from the SEEDED fits
_lq = np.log(QV)
_sl = np.array([np.polyfit(_lq, np.log(r[4]), 1)[0] for r in S['seeded']
                if np.all(np.isfinite(r[4])) and np.all(np.asarray(r[4]) > 0)])
_f0 = np.array([fd[2][Q0] for fd, r in zip(S['diff'], S['seeded'])
                if np.all(np.isfinite(r[4])) and np.all(np.asarray(r[4]) > 0)])
_k = _f0 > 0.15
print(f'\n@@ Q exponent of the released tau_fast (seeded): median {np.median(_sl[_k]):.2f} over the '
      f'{_k.sum()} groups with f > 0.15 at the lowest Q, 80% between {np.percentile(_sl[_k], 10):.1f} '
      f'and {np.percentile(_sl[_k], 90):.1f}')
print(f'\n@@ groups where the DEFAULT seed lands above chi2 10 ({len(bad)}):')
for run, tw, cd, cs in bad:
    print(f'@@   {run} t_w={tw:.0f}s  default chi2 {cd:.1f}  vs seeded {cs:.2f}')

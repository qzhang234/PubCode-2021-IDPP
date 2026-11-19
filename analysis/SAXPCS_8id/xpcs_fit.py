"""The two-mode correlation-function model (Equation 2 of the main text) and its fits.

Every figure script and general_model_check.py fit g2 through the functions
here, so all figures show one model applied one way.  The choices behind the
fit are described in the repository README.
"""

import h5py
import numpy as np
from scipy.optimize import least_squares

from nexus_read import read_g2

# Speckle contrast beta.  The calling script measures it on the static standard
# and sets it before fitting:
#     xpcs_fit.CONTRAST = xpcs_fit.measure_beta(STANDARD, label)
CONTRAST = None
BASELINE = 1.0

# beta is read from the flat part of the static standard, not from the five
# fitted bins: the 2021-3 standard (NIST glassy carbon) is scratched, and its
# static parasitic scattering depresses the apparent contrast at low Q, by 12%
# in the lowest fitted bin.  The band stops before a dip at 0.0175-0.0185 A^-1
# and before the bin of the stuck pixel at 0.0214 A^-1.
BETA_BAND = (0.0085, 0.0165)          # A^-1

# p_fast values of the profile that seeds every fit (fit_seeded)
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)

# Start values and bounds.  Per Q bin: (tau_fast [s], f, tau_slow [s]).
# Shared by the bins of one group: the stretching exponents (p_fast, p_slow).
PQ_P0    = [1e-3, 0.5, 100.0]
PQ_LO    = [1e-6, 0.0, 1.0]
PQ_HI    = [10.0, 1.0, 10000.0]
P_EXP_P0 = [0.5, 0.5]
P_EXP_LO = [0.2, 0.2]
P_EXP_HI = [3.0, 3.0]

# D0 is of order 1e7 A^2/s while every other parameter is of order 1, so it is
# fitted in units of D_SCALE; otherwise its column of the Jacobian is so small
# that the pseudo-inverse returns a zero uncertainty for it.
D_SCALE = 1e7
D_LO, D_HI = 1e-5, 1e3          # in units of D_SCALE


def flat_level(tau, g, e):
    """Inverse-variance weighted mean of g2 over the measured delays, and its
    standard error."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * g[ok]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


def measure_beta(path, label, fitted_bins=(0, 1, 2, 3, 4)):
    """beta of a static reference: g2 of a static sample is flat at 1 + beta.

    The flat level is taken in every Q bin, and beta is the inverse-variance
    mean of level - 1 over the bins inside BETA_BAND.
    """
    with h5py.File(path, 'r') as hf:
        tau, g2, g2e, q_std = read_g2(hf)
    n_bins = g2.shape[1]
    b = np.zeros(n_bins)
    e = np.zeros(n_bins)
    for i in range(n_bins):
        level, level_err = flat_level(tau, g2[:, i], g2e[:, i])
        b[i] = level - 1.0
        e[i] = level_err
    sel = (q_std >= BETA_BAND[0]) & (q_std <= BETA_BAND[1]) & np.isfinite(b)
    w = 1.0 / e[sel] ** 2
    beta = float(np.sum(w * b[sel]) / np.sum(w))
    print(f'beta ({label}) = {beta:.5f} from {sel.sum()} bins over '
          f'{BETA_BAND[0]}-{BETA_BAND[1]} A^-1 (spread {b[sel].min():.4f}-{b[sel].max():.4f})')
    fitted_text = []
    for i in fitted_bins:
        fitted_text.append(f'{b[i]:.4f}')
    print('   over the five FITTED bins it would read ' + ', '.join(fitted_text)
          + '  <- depressed by the standard')
    return beta


def double_exp(tau, tau_fast, f, tau_slow, p1, p2):
    """Equation 2: g2 = 1 + beta * [f g1_fast + (1 - f) g1_slow]^2.

    Each mode is a stretched exponential, exp[-(tau/tau_i)^p_i].  The modes are
    added as field correlations and the sum is squared (Siegert relation).
    """
    if CONTRAST is None:
        raise ValueError('xpcs_fit.CONTRAST is not set: measure beta with measure_beta() first')
    decay_fast = f * np.exp(-(tau / tau_fast)**p1)
    decay_slow = (1 - f) * np.exp(-(tau / tau_slow)**p2)
    return CONTRAST * (decay_fast + decay_slow)**2 + BASELINE


def usable_bins(tau, g2, g2_err, q_indices):
    """(bin, tau, g2, g2_err) for every bin in q_indices with at least five
    usable delays: positive delay, finite g2 and an uncertainty above zero."""
    data = []
    for qi in q_indices:
        usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
        if usable.sum() >= 5:
            data.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
    return data


def bin_params(p, i, offset):
    """(tau_fast, f, tau_slow) of the i-th Q bin in the parameter vector p.

    The first `offset` entries of p are the shared exponents; after them come
    three numbers per bin, so bin i starts at position offset + 3 * i.
    """
    start = offset + 3 * i
    return p[start], p[start + 1], p[start + 2]


def covariance_errors(res):
    """One-standard-deviation errors from the least-squares Jacobian.

    The covariance is the pseudo-inverse of J^T J, so a parameter that settles
    on a bound (f = 1 happens in some mobile groups) leaves the uncertainties
    of the other parameters defined.
    """
    cov = np.linalg.pinv(res.jac.T @ res.jac, rcond=1e-12)
    return np.sqrt(np.abs(np.diag(cov)))


def fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=None, p1_start=P_EXP_P0[0]):
    """Two-mode fit of several Q bins at once, with tau_fast free in every bin.

    This is the general two-mode form; Equation 2 is fit_g2_diffusive().  The
    exponents p_fast (p1) and p_slow (p2) are shared by the bins; p1_fixed
    holds p_fast at a given value, and p1_start is where a free p_fast starts.
    Minimises the error-weighted residuals (model - g2) / g2_err.

    Returns {'p1', 'p1_err', 'p2', 'p2_err', 'red_chi2', 'n_par',
    'per_q': {bin: {'tau_fast', 'tau_fast_err', 'f', 'f_err', 'tau_slow',
    'tau_slow_err'}}}, or None if no bin has usable data.
    """
    data = usable_bins(tau, g2, g2_err, q_indices)
    nq = len(data)
    if nq == 0:
        return None

    # parameter vector: [p_fast (only when free), p_slow] + (tau_fast, f, tau_slow) per bin
    if p1_fixed is None:
        offset = 2
        x0 = [p1_start, P_EXP_P0[1]]
        lo = [P_EXP_LO[0], P_EXP_LO[1]]
        hi = [P_EXP_HI[0], P_EXP_HI[1]]
    else:
        offset = 1
        x0 = [P_EXP_P0[1]]
        lo = [P_EXP_LO[1]]
        hi = [P_EXP_HI[1]]
    for i in range(nq):
        x0 = x0 + PQ_P0
        lo = lo + PQ_LO
        hi = hi + PQ_HI

    def residual(p):
        if p1_fixed is None:
            p1, p2 = p[0], p[1]
        else:
            p1, p2 = p1_fixed, p[0]
        parts = []
        for i in range(nq):
            qi, tv, gv, ev = data[i]
            tf, f, ts = bin_params(p, i, offset)
            parts.append((double_exp(tv, tf, f, ts, p1, p2) - gv) / ev)
        return np.concatenate(parts)

    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=40000)
    ndof = max(len(res.fun) - len(res.x), 1)
    perr = covariance_errors(res)

    if p1_fixed is None:
        p1v, p1e = res.x[0], perr[0]
        p2v, p2e = res.x[1], perr[1]
    else:
        p1v, p1e = p1_fixed, 0.0
        p2v, p2e = res.x[0], perr[0]
    out = {'p1': p1v, 'p1_err': p1e, 'p2': p2v, 'p2_err': p2e,
           'red_chi2': float(np.sum(res.fun**2) / ndof), 'n_par': len(res.x), 'per_q': {}}
    for i in range(nq):
        qi = data[i][0]
        tf, f, ts = bin_params(res.x, i, offset)
        tf_err, f_err, ts_err = bin_params(perr, i, offset)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf_err,
                            'f': f, 'f_err': f_err,
                            'tau_slow': ts, 'tau_slow_err': ts_err}
    return out


def fit_g2_diffusive(tau, g2, g2_err, q_indices, q_vals, p1_fixed=None, p1_start=P_EXP_P0[0]):
    """Equation 2: the two-mode fit with tau_fast(Q) = 1 / (D0 Q^2).

    The Q^-2 form is imposed, not fitted.  At delay 1/(D0 Q^2) the fast term
    has fallen to 1/e of its amplitude whatever p_fast is, so D0 is the rate of
    the fast decay.  p1_fixed and p1_start are as in fit_g2_global(), which
    supplies the starting values.

    Returns a fit_g2_global-style dict with 'D0' and 'D0_err' (A^2/s) added,
    or None if no bin has usable data.
    """
    data = usable_bins(tau, g2, g2_err, q_indices)
    nq = len(data)
    if nq == 0:
        return None
    qs = np.zeros(nq)
    for i in range(nq):
        qs[i] = q_vals[data[i][0]]

    # start from the general fit; its tau_fast values give a starting D0
    free = fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=p1_fixed, p1_start=p1_start)
    if free is None:
        return None
    d_values = []
    for i in range(nq):
        qi = data[i][0]
        d_values.append(1.0 / (free['per_q'][qi]['tau_fast'] * qs[i] ** 2))
    d_seed = float(np.median(d_values)) / D_SCALE

    # parameter vector: [p_fast (only when free), p_slow, D0 / D_SCALE] + (f, tau_slow) per bin
    if p1_fixed is None:
        head = 3
        x0 = [free['p1'], free['p2'], d_seed]
        lo = [P_EXP_LO[0], P_EXP_LO[1], D_LO]
        hi = [P_EXP_HI[0], P_EXP_HI[1], D_HI]
    else:
        head = 2
        x0 = [free['p2'], d_seed]
        lo = [P_EXP_LO[1], D_LO]
        hi = [P_EXP_HI[1], D_HI]
    i_pslow = head - 2
    i_d0 = head - 1
    for i in range(nq):
        pq = free['per_q'][data[i][0]]
        x0 = x0 + [pq['f'], pq['tau_slow']]
        lo = lo + [PQ_LO[1], PQ_LO[2]]
        hi = hi + [PQ_HI[1], PQ_HI[2]]
    x0 = np.clip(x0, lo, hi)

    def residual(p):
        if p1_fixed is None:
            p_fast = p[0]
        else:
            p_fast = p1_fixed
        p_slow = p[i_pslow]
        D0 = p[i_d0] * D_SCALE
        parts = []
        for i in range(nq):
            qi, tv, gv, ev = data[i]
            f = p[head + 2 * i]
            ts = p[head + 2 * i + 1]
            parts.append((double_exp(tv, 1.0 / (D0 * qs[i] ** 2), f, ts, p_fast, p_slow)
                          - gv) / ev)
        return np.concatenate(parts)

    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=120000)
    ndof = max(len(res.fun) - len(res.x), 1)
    perr = covariance_errors(res)
    if p1_fixed is None:
        p_fast, p_fast_err = res.x[0], perr[0]
    else:
        p_fast, p_fast_err = p1_fixed, 0.0
    D0 = res.x[i_d0] * D_SCALE
    D0_err = perr[i_d0] * D_SCALE
    out = {'p1': p_fast, 'p1_err': p_fast_err,
           'p2': res.x[i_pslow], 'p2_err': perr[i_pslow],
           'D0': D0, 'D0_err': D0_err,
           'red_chi2': float(np.sum(res.fun**2) / ndof),
           'n_par': len(res.x), 'per_q': {}}
    for i in range(nq):
        qi = data[i][0]
        tf = 1.0 / (D0 * qs[i] ** 2)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf * D0_err / D0,
                            'f': res.x[head + 2 * i], 'f_err': perr[head + 2 * i],
                            'tau_slow': res.x[head + 2 * i + 1],
                            'tau_slow_err': perr[head + 2 * i + 1]}
    return out


def fit_model(model, tau, g2, g2_err, q_indices, q_vals, p1_fixed=None, p1_start=P_EXP_P0[0]):
    """Equation 2 (model='diffusive') or the general two-mode form (model='general')."""
    if model == 'diffusive':
        return fit_g2_diffusive(tau, g2, g2_err, q_indices, q_vals,
                                p1_fixed=p1_fixed, p1_start=p1_start)
    if model == 'general':
        return fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=p1_fixed, p1_start=p1_start)
    raise ValueError(f'unknown model {model!r}')


def fit_seeded(tau, g2, g2_err, q_indices, q_vals, model='diffusive'):
    """The fit behind every figure: p_fast free, started from a p_fast profile.

    With p_fast free the two-mode fit has local minima in which the fast term
    imitates the slow one.  So the fit is first run with p_fast held at each
    value of P_GRID, then released from the best of those, and the released fit
    is kept only if its chi^2 is no worse.  Returns None if no bin has usable data.
    """
    best = None
    chi_best = None
    p_best = None
    for p_fast in P_GRID:
        r = fit_model(model, tau, g2, g2_err, q_indices, q_vals, p1_fixed=p_fast)
        if r is None:
            continue
        if best is None or r['red_chi2'] < chi_best:
            best, chi_best, p_best = r, r['red_chi2'], p_fast
    if best is None:
        return None
    free = fit_model(model, tau, g2, g2_err, q_indices, q_vals, p1_start=p_best)
    if free is not None and free['red_chi2'] <= chi_best:
        return free
    return best

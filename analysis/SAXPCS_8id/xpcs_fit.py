"""Shared XPCS two-mode fit: the model, the measured contrast, and the global fit.

Figure 3b,c, Figure S9 and Figure S10 must all show the SAME model applied the
SAME way -- Figure S10 is an expansion of Figure 3b onto the four higher q bins
-- so the model and the fitting routine live here and every figure script
imports them.  (Before this module, g2_grid_SI.py carried its own copy that
fitted each q bin independently with the stretching exponents frozen at 0.5, so
Figure S10 was showing a different model from the one its caption described.)

CONTRAST
--------
beta is the instrumental speckle contrast: a property of the beam coherence and
the detector pixel size, not of the sample.  It is therefore measured once, on a
static reference, and held fixed for every sample fit.

The value below was measured on the 10 nm nano-porous glass standard (Doraglas
S10-10-1200-50, 10 mm x 1.2 mm, 10 nm pores), dataset F0145, at 6 C: 50 repeat
acquisitions, q bin 16 (q = 0.02067 A^-1, the strongest bin whose averaged
correlation function is flat).  Averaging the 50 repeats and fitting with a straight
line in log(delay) gives a flat correlation function -- slope -3.1(1.8)e-5, i.e.
1.7 sigma from zero, reduced chi^2 = 1.02 -- as a static sample must, with

    g2 = 1 + beta,   beta = 0.13042 +/- 0.00001.

All 50 repeats were retained: at this q bin every one is statistically flat
(worst reduced chi^2 against a constant = 1.30).  The visibly non-flat repeats
are at the neighbouring bin 15 (q = 0.01954 A^-1), where the median reduced
chi^2 against a constant is 12.3; that bin is not used.
"""

import numpy as np
from scipy.optimize import least_squares

# --- measured instrumental contrast (see the module docstring) ---
CONTRAST = 0.13042
BASELINE = 1.0

# per-q parameter bounds / start (tau_fast, f, tau_slow) and shared (p1, p2)
PQ_P0    = [1e-3, 0.5, 100.0]
PQ_LO    = [1e-6, 0.0, 1.0]
PQ_HI    = [10.0, 1.0, 10000.0]
P_EXP_P0 = [0.5, 0.5]
P_EXP_LO = [0.2, 0.2]
P_EXP_HI = [3.0, 3.0]


def double_exp(tau, tau_fast, f, tau_slow, p1, p2):
    decay_fast = f * np.exp(-(tau / tau_fast)**p1)
    decay_slow = (1 - f) * np.exp(-(tau / tau_slow)**p2)
    return CONTRAST * (decay_fast + decay_slow)**2 + BASELINE


def fit_g2_global(tau, g2, g2_err, q_indices):
    """Global fit of several q bins for one elapsed time, sharing p1 and p2.

    Parameter vector = [p1, p2, (tau_fast, f, tau_slow) x nq].  Minimises the
    error-weighted residual (model - g2) / g2_err over all q simultaneously,
    using the g2_err stored in the file directly (absolute_sigma convention).
    Parameter 1-sigma errors are sqrt(diag(inv(J^T J))).

    Returns a dict:
      {'p1','p1_err','p2','p2_err','red_chi2',
       'per_q': {q_idx: {'tau_fast','tau_fast_err','f','f_err',
                         'tau_slow','tau_slow_err'}}}
    or None if too few q bins have usable data.
    """
    data = []
    for qi in q_indices:
        v = (tau > 0) & ~np.isnan(g2[:, qi]) & ~np.isnan(g2_err[:, qi]) & (g2_err[:, qi] > 0)
        if v.sum() >= 5:
            data.append((qi, tau[v], g2[v, qi], g2_err[v, qi]))
    nq = len(data)
    if nq == 0:
        return None

    def residual(p):
        p1, p2 = p[0], p[1]
        parts = []
        for i, (qi, tv, gv, ev) in enumerate(data):
            tf, f, ts = p[2 + 3 * i: 5 + 3 * i]
            parts.append((double_exp(tv, tf, f, ts, p1, p2) - gv) / ev)
        return np.concatenate(parts)

    x0 = list(P_EXP_P0) + PQ_P0 * nq
    lo = list(P_EXP_LO) + PQ_LO * nq
    hi = list(P_EXP_HI) + PQ_HI * nq
    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=40000)

    ndof = max(len(res.fun) - len(res.x), 1)
    red_chi2 = float(np.sum(res.fun**2) / ndof)
    # covariance from the Gauss-Newton Hessian of error-weighted residuals
    # Pseudo-inverse, not inverse: at the latest elapsed time f rails against its
    # lower bound at the lowest q, which makes tau_fast unidentifiable there and
    # J^T J exactly singular.  inv() then returns nan for EVERY parameter,
    # including the well-constrained ones (p1, p2 and the other q bins).  pinv
    # discards only the degenerate direction, so the identifiable parameters keep
    # honest uncertainties and only the unidentifiable one is reported as ~0.
    jtj = res.jac.T @ res.jac
    cov = np.linalg.pinv(jtj, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))

    out = {'p1': res.x[0], 'p1_err': perr[0],
           'p2': res.x[1], 'p2_err': perr[1],
           'red_chi2': red_chi2, 'per_q': {}}
    for i, (qi, tv, gv, ev) in enumerate(data):
        tf, f, ts = res.x[2 + 3 * i: 5 + 3 * i]
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': perr[2 + 3 * i],
                            'f': f, 'f_err': perr[3 + 3 * i],
                            'tau_slow': ts, 'tau_slow_err': perr[4 + 3 * i]}
    return out

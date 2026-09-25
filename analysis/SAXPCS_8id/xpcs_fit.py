"""Shared XPCS two-mode fit: the model, the measured contrast, and the global fit.

Figure 3b,c, Figure S8 and Figure S9 must all show the SAME model applied the
SAME way -- Figure S9 is an expansion of Figure 3b onto the four higher q bins
-- so the model and the fitting routine live here and every figure script
imports them.  (Before this module, g2_grid_SI.py carried its own copy that
fitted each q bin independently with the stretching exponents frozen at 0.5, so
Figure S9 was showing a different model from the one its caption described.)

CONTRAST
--------
beta is the instrumental speckle contrast: a property of the beam coherence and
the detector pixel size, not of the sample.  It is therefore measured once, on a
static reference, and held fixed for every sample fit.

The value below was measured on the 10 nm nano-porous glass standard (Doraglas
S10-10-1200-50, 10 mm x 1.2 mm, 10 nm pores), dataset F0145, at 6 C.  The value
comes from the group average average_ranges.py writes into data/ -- 41 of the 50
acquisitions, the survivors of the same two outlier cuts every other group in
the paper gets -- and it is measured in the SAME five q bins the sample is
fitted in, which is where it is used.  A constant is fitted to the averaged g2
of each, as a static sample requires; the five are flat (reduced chi^2 0.71 to
1.03) and return 0.12948, 0.13124, 0.13291, 0.13145 and 0.13186, whose mean is

    g2 = 1 + beta,   beta = 0.13139 +/- 0.00056,

the error being the standard error of that mean.  The strongest-counted flat bin
anywhere, bin 16 at q = 0.02067 A^-1, gives 0.13031, within one standard error,
so nothing hinges on the choice: moving between the two shifts every fitted fast
fraction by less than 0.02.
"""

import numpy as np
from scipy.optimize import least_squares

# --- measured instrumental contrast (see the module docstring) ---
CONTRAST = 0.13139
BASELINE = 1.0

# per-q parameter bounds / start (tau_fast, f, tau_slow) and shared (p1, p2)
PQ_P0    = [1e-3, 0.5, 100.0]
PQ_LO    = [1e-6, 0.0, 1.0]
PQ_HI    = [10.0, 1.0, 10000.0]
P_EXP_P0 = [0.5, 0.5]
P_EXP_LO = [0.2, 0.2]
P_EXP_HI = [3.0, 3.0]


def double_exp(tau, tau_fast, f, tau_slow, p1, p2):
    """Two-mode g2: a fraction f of the scattering relaxes fast, 1-f slowly.

    Each mode is a stretched exponential (Kohlrausch-Williams-Watts): p = 1 is
    a plain exponential, p < 1 stretches the decay, p > 1 compresses it.  The
    two are added as FIELDS and then squared, because g2 is the intensity
    correlation and the field correlation is what decays (Siegert relation).
    """
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
    or None if no q bin has usable data.
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
            # p = [p1, p2, then three numbers per q bin], so bin i owns
            # p[2 + 3i], p[3 + 3i], p[4 + 3i].
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

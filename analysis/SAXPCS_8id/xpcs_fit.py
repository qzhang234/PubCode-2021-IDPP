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


def bin_params(p, i):
    """The three parameters belonging to q bin number i.

    The fit works on one long list of numbers, p.  The first two entries are
    the shared exponents p1 and p2.  After those come three numbers for each q
    bin, always in the order tau_fast, f, tau_slow.  So bin 0 starts at
    position 2, bin 1 at position 5, and bin i at position 2 + 3 * i.
    """
    start = 2 + 3 * i
    return p[start], p[start + 1], p[start + 2]


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
        # A delay point is usable when the delay is positive, g2 is a real
        # number, and its uncertainty is a real number above zero, since the
        # fit divides by that uncertainty.
        usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
        if usable.sum() >= 5:
            data.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
    nq = len(data)
    if nq == 0:
        return None

    def residual(p):
        p1, p2 = p[0], p[1]
        parts = []
        for i, (qi, tv, gv, ev) in enumerate(data):
            tf, f, ts = bin_params(p, i)
            parts.append((double_exp(tv, tf, f, ts, p1, p2) - gv) / ev)
        return np.concatenate(parts)

    # "PQ_P0 * nq" repeats that list of three numbers once per q bin, which
    # builds the starting guess for the whole parameter vector in one line.
    x0 = list(P_EXP_P0) + PQ_P0 * nq
    lo = list(P_EXP_LO) + PQ_LO * nq
    hi = list(P_EXP_HI) + PQ_HI * nq
    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=40000)

    ndof = max(len(res.fun) - len(res.x), 1)
    red_chi2 = float(np.sum(res.fun**2) / ndof)
    # Parameter uncertainties from the Gauss-Newton approximation to the
    # Hessian of the error-weighted residuals.
    #
    # Pseudo-inverse rather than plain inverse, as a safeguard.  If an amplitude
    # f ever settles exactly on its bound of 0, its tau_fast stops affecting the
    # model, J^T J becomes singular, and inv() returns nan for EVERY parameter,
    # including the well-determined ones.  pinv drops only the useless direction
    # and leaves the rest with honest uncertainties.  No parameter reaches a
    # bound at the contrast used here, so the two agree in practice.
    jtj = res.jac.T @ res.jac
    cov = np.linalg.pinv(jtj, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))

    out = {'p1': res.x[0], 'p1_err': perr[0],
           'p2': res.x[1], 'p2_err': perr[1],
           'red_chi2': red_chi2, 'per_q': {}}
    for i, (qi, tv, gv, ev) in enumerate(data):
        tf, f, ts = bin_params(res.x, i)
        tf_err, f_err, ts_err = bin_params(perr, i)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf_err,
                            'f': f, 'f_err': f_err,
                            'tau_slow': ts, 'tau_slow_err': ts_err}
    return out

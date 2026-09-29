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
from scipy.special import gamma as gamma_fn

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


def bin_params(p, i, offset=2):
    """The three parameters belonging to q bin number i.

    The fit works on one long list of numbers, p.  The first `offset` entries
    are shared: the exponents p1 and p2 when they are free, nothing when they
    are held fixed.  After those come three numbers for each q bin, always in
    the order tau_fast, f, tau_slow.  So with the default offset bin 0 starts
    at position 2, bin 1 at position 5, and bin i at position 2 + 3 * i.
    """
    start = offset + 3 * i
    return p[start], p[start + 1], p[start + 2]


def fit_g2_global(tau, g2, g2_err, q_indices, p1_fixed=None):
    """Global fit of several q bins for one elapsed time, sharing p1 and p2.

With p1_fixed the FAST stretching exponent is held at the value fit_g2_joint
    returns for the whole series, so that tau_fast is comparable between elapsed
    times; the parameter vector then drops its first entry.  p_slow stays free
    at every elapsed time -- tau_slow lies beyond the acquisition window and is
    not used for scaling, so there is nothing to gain by constraining it and the
    fit is better without.

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

    # Fixing p_fast removes one entry from the head of the parameter vector, so
    # every per-q block shifts down by one and the offset is 1 instead of 2.
    off = 1 if p1_fixed is not None else 2

    def residual(p):
        p1, p2 = (p1_fixed, p[0]) if p1_fixed is not None else (p[0], p[1])
        parts = []
        for i, (qi, tv, gv, ev) in enumerate(data):
            tf, f, ts = bin_params(p, i, offset=off)
            parts.append((double_exp(tv, tf, f, ts, p1, p2) - gv) / ev)
        return np.concatenate(parts)

    # "PQ_P0 * nq" repeats that list of three numbers once per q bin, which
    # builds the starting guess for the whole parameter vector in one line.
    # head = [p_slow] when p_fast is fixed, [p_fast, p_slow] when neither is
    x0 = (P_EXP_P0[1:] if p1_fixed is not None else list(P_EXP_P0)) + PQ_P0 * nq
    lo = (P_EXP_LO[1:] if p1_fixed is not None else list(P_EXP_LO)) + PQ_LO * nq
    hi = (P_EXP_HI[1:] if p1_fixed is not None else list(P_EXP_HI)) + PQ_HI * nq
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

    if p1_fixed is not None:
        p1v, p1e = p1_fixed, 0.0
        p2v, p2e = res.x[0], perr[0]
    else:
        p1v, p1e = res.x[0], perr[0]
        p2v, p2e = res.x[1], perr[1]
    out = {'p1': p1v, 'p1_err': p1e, 'p2': p2v, 'p2_err': p2e,
           'red_chi2': red_chi2, 'n_par': len(res.x), 'per_q': {}}
    for i, (qi, tv, gv, ev) in enumerate(data):
        tf, f, ts = bin_params(res.x, i, offset=off)
        tf_err, f_err, ts_err = bin_params(perr, i, offset=off)
        out['per_q'][qi] = {'tau_fast': tf, 'tau_fast_err': tf_err,
                            'f': f, 'f_err': f_err,
                            'tau_slow': ts, 'tau_slow_err': ts_err}
    return out


def mean_tau(tau_se, p):
    """Mean relaxation time of a stretched exponential.

    A Kohlrausch decay exp[-(t/tau_se)^p] does not relax on the time tau_se;
    its mean is

        <tau> = tau_se * Gamma(1/p) / p,

    which is the quantity to compare when p is not the same everywhere (Guo et
    al., Phys. Rev. Lett. 109, 055901 (2012), Eq. 1 and the text below it).
    With p held global across the series this is a constant rescaling, so the
    q dependence is unchanged and only the absolute time shifts; the reason to
    apply it is that <tau> is the physically meaningful relaxation time and is
    what the XPCS literature quotes.
    """
    return tau_se * gamma_fn(1.0 / p) / p


def fit_g2_joint(datasets, q_indices):
    """One fit of the whole waiting-time series, sharing p_fast across it.

    WHY.  Fitted separately, each elapsed time returns its own p_fast, and those
    ran 0.47 to 0.68 across the series.  A Kohlrausch tau is not comparable
    between fits with different p -- the mean time carries a factor
    Gamma(1/p)/p that moves by 50 % over that range -- so a trend in tau_fast
    against waiting time could not be read as a trend in relaxation rate.
    p_fast is therefore fitted once for the series.  One value describes it: the
    chi^2/dof profile is flat to 1 % between 0.55 and 0.70, and letting it float
    per time improves chi^2/dof only from 1.201 to 1.188 for four more
    parameters.

    p_slow is left free at every elapsed time.  tau_slow lies beyond the
    acquisition window and is never used for scaling, so constraining it buys
    nothing and costs fit quality.

    tau_fast stays free in every q bin, so the diffusive q^-2 scaling remains a
    result to be tested rather than an assumption; it holds within uncertainty
    at the elapsed times where the fast mode is resolved.

    datasets: list of (tau, g2, g2_err) in waiting-time order.
    Returns (joint, per_time):
      joint    = {'p1','p1_err','red_chi2','n_par'}      # p1 = shared p_fast
      per_time = list of fit_g2_global-style dicts, one per dataset, refitted at
                 the shared p_fast so callers get the familiar per_q structure
                 and each time keeps its own p_slow.
    """
    # Seed from independent per-time fits.
    seeds = [fit_g2_global(t, g, e, q_indices) for t, g, e in datasets]
    if not any(s is not None for s in seeds):
        return None, None
    p1_0 = float(np.mean([s['p1'] for s in seeds if s is not None]))

    # Assemble the usable delay points once, in the order the residual walks them.
    blocks = []
    for (tau, g2, g2_err) in datasets:
        d = []
        for qi in q_indices:
            usable = (tau > 0) & np.isfinite(g2[:, qi]) & (g2_err[:, qi] > 0)
            if usable.sum() >= 5:
                d.append((qi, tau[usable], g2[usable, qi], g2_err[usable, qi]))
        blocks.append(d)

    # Parameter vector: [p_fast] then, per elapsed time, [p_slow, (tau_fast, f,
    # tau_slow) x nq].  starts[k] is where elapsed time k's p_slow sits, so its
    # q bins begin one entry later.
    starts, x0, lo, hi = [], [p1_0], [0.2], [3.0]
    for d, s in zip(blocks, seeds):
        starts.append(len(x0))
        x0.append(s['p2'] if s else P_EXP_P0[1])
        lo.append(P_EXP_LO[1])
        hi.append(P_EXP_HI[1])
        for i, (qi, _, _, _) in enumerate(d):
            pq = s['per_q'].get(qi) if s else None
            x0 += [pq['tau_fast'], pq['f'], pq['tau_slow']] if pq else list(PQ_P0)
            lo += list(PQ_LO)
            hi += list(PQ_HI)
    x0 = [min(max(v, l), h) for v, l, h in zip(x0, lo, hi)]

    def residual(p):
        parts = []
        for d, st in zip(blocks, starts):
            p2 = p[st]
            for i, (qi, tv, gv, ev) in enumerate(d):
                tf, f, ts = bin_params(p, i, offset=st + 1)
                parts.append((double_exp(tv, tf, f, ts, p[0], p2) - gv) / ev)
        return np.concatenate(parts)

    res = least_squares(residual, x0, bounds=(lo, hi), max_nfev=300000)
    ndof = max(len(res.fun) - len(res.x), 1)
    cov = np.linalg.pinv(res.jac.T @ res.jac, rcond=1e-12)
    perr = np.sqrt(np.abs(np.diag(cov)))
    joint = {'p1': res.x[0], 'p1_err': perr[0],
             'red_chi2': float(np.sum(res.fun**2) / ndof), 'n_par': len(res.x)}

    # Refit each elapsed time at the shared p_fast.  Same parameters the joint
    # fit found, but each time then carries its own covariance and p_slow, and
    # callers get the per_q dict they already expect.
    per_time = [fit_g2_global(t, g, e, q_indices, p1_fixed=res.x[0])
                for t, g, e in datasets]
    return joint, per_time

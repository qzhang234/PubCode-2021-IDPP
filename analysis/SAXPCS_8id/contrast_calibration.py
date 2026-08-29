"""Figure S4: measurement of the instrumental speckle contrast beta.

beta is a property of the beamline -- the transverse coherence of the incident
beam and the speckle-to-pixel size ratio -- not of the sample, so it is measured
once on a STATIC reference and then held fixed in every g2 fit.  A static sample
has no decorrelation, so its correlation function is flat at

    g2(tau) = 1 + beta.

Reference: nano-porous glass, 10 nm nominal pore diameter, 10 x 10 mm, 1.2 mm
thick (Doraglas S10-10-1200-50).  Dataset F0145, 6 C, 50 repeat acquisitions of
100,000 frames each, taken at 19:49 on 2022-03-06 -- about two hours before the
isothermal series of Figure 3 (21:57), in the same beamline configuration.  The
50 correlation functions are read from the stack average_ranges.py writes into
data/, so this script does not touch the beamline storage either.

Panel (a): the 50 repeats and their average at q = 0.02067 A^-1, with the
straight-line fit whose level gives beta.  That bin is not the best-counted one.
Ranking all 27 bins by mean pointwise g2_err -- the committed stack carries no
intensity dataset, so g2_err is the available proxy for counting statistics --
bin 16 ranks 8th, and every one of the seven bins ahead of it departs further
from flat: reduced chi^2 of 1.87, 2.20, 2.67, 9.61, 3.86, 2.43 and 17.09 against
a constant, versus 1.02 here.  Bin 16 is therefore the best-counted bin whose
averaged g2 is flat to within its own uncertainties, and so measures the
instrument rather than residual structure in the standard.  Panel (b): beta from
the same procedure at every q bin, showing that it is essentially independent of
q, so the value measured there transfers to the low-q bins used for the sample.
Four of those five bins also carry their own direct measurement: beta = 0.1297,
0.1328, 0.1318 and 0.1319 at q = 0.00376, 0.00602, 0.00714 and 0.00827 A^-1, all
within 2 % of the fixed 0.13042.

Four of the 27 bins are excluded from panel (b) and flagged in the print-out
(bins 2, 12, 15 and 22): their averaged correlation function is not flat,
reduced chi^2 > 5 against a constant.  Bin 2 (q = 0.00489 A^-1) is the extreme
case, returning beta = 0.65 with a reduced chi^2 of 2346.  These are
detector/parasitic-scattering artefacts, not contrast.  Every one of the 50
repeats at the bin actually used is statistically flat (worst reduced chi^2 =
1.30), so no repeat is rejected.

Panel (b) is flat at about 0.132 out to q ~ 0.021 A^-1 and then falls off,
reaching 0.106 in the outermost bin at q = 0.0331 A^-1.  That falloff is real
and instrumental, not a defect of the standard: the 0.03 % bandwidth passed by
the Si(111) monochromator at 10.91 keV gives a finite LONGITUDINAL coherence
length, and the path-length difference across the scattering volume grows with
scattering angle, so contrast is progressively lost as q rises.  It is the
reason bin 16 (q = 0.02067 A^-1) is a good choice and a bin near the top of the
range would not be: bin 16 still sits on the plateau, and so do all five sample
bins, which are at q <= 0.00827 A^-1.
"""

import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)

# The 50 repeats, stacked into one file by average_ranges.py.  They are stacked
# rather than averaged because panel (a) plots every one of them: the scatter
# between repeats is part of what the figure shows.  The stack keeps the raw
# NeXus paths, with the repeat index as the leading axis of g2 and g2_err.
GLASS = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data',
                     'Stack_F0145_10nm_Glass_006C_att00_Rq0_00001_00050_results.hdf')
Q_INDEX = 15                 # 0-based (bin 16); q = 0.02067 A^-1, the best-counted
                             # bin whose averaged g2 is flat (chi^2_red = 1.02)
CHI2_FLAT_MAX = 5.0          # a q bin whose averaged g2 is this far from flat is
                             # an artefact, not a contrast measurement
FIG_SIZE = (DOUBLE_COL, 2.9)


def load(path):
    """Delay times, q bins, and the (repeat, delay, q) g2 and g2_err stacks."""
    with h5py.File(path, 'r') as hf:
        ft = float(np.asarray(hf['/entry/instrument/detector_1/frame_time'][()]).reshape(-1)[0])
        t = np.asarray(hf['/xpcs/multitau/delay_list'][()])
        tau = (t[:, 0] if t.ndim > 1 else t) * ft
        q = hf['/xpcs/qmap/dynamic_v_list_dim0'][()]
        G = hf['/xpcs/multitau/normalized_g2'][()]
        E = hf['/xpcs/multitau/normalized_g2_err'][()]
        n = hf['/xpcs/average/file_list'].shape[0]
    return tau, np.asarray(q), np.asarray(G), np.asarray(E), n


def flat_level(tau, g, e):
    """Weighted mean level of a flat g2, its error, and the reduced chi^2 of the
    flat hypothesis."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & np.isfinite(e) & (e > 0)
    w = 1.0 / e[ok]**2
    lvl = np.sum(w * g[ok]) / np.sum(w)
    return lvl, 1.0 / np.sqrt(np.sum(w)), np.sum(w * (g[ok] - lvl)**2) / (ok.sum() - 1)


tau, q, G, E, n_rep = load(GLASS)
print(f'{n_rep} repeats of the 10 nm glass standard, {len(q)} q bins')

# --- per-repeat check at the working bin: is any repeat an outlier? ---
worst = max(flat_level(tau, G[i, :, Q_INDEX], E[i, :, Q_INDEX])[2] for i in range(n_rep))
print(f'q bin {Q_INDEX+1} (q = {q[Q_INDEX]:.5f} A^-1): worst per-repeat '
      f'reduced chi^2 against a constant = {worst:.2f} -> no repeat rejected')

# --- averaged g2 at the working bin, and the straight-line fit ---
gbar = np.nanmean(G[:, :, Q_INDEX], axis=0)
ebar = np.sqrt(np.nansum(E[:, :, Q_INDEX]**2, axis=0)) / n_rep
ok = (tau > 0) & (tau < 2) & np.isfinite(gbar) & (ebar > 0)
x, y, w = np.log10(tau[ok]), gbar[ok], 1.0 / ebar[ok]**2
X = np.vstack([x, np.ones_like(x)]).T
cov = np.linalg.inv(X.T @ (w[:, None] * X))
slope, intercept = cov @ (X.T @ (w * y))
beta, beta_err, chi2 = flat_level(tau, gbar, ebar)
beta -= 1.0
print(f'  straight-line slope = {slope:+.2e} +/- {np.sqrt(cov[0,0]):.2e} '
      f'({abs(slope)/np.sqrt(cov[0,0]):.1f} sigma from flat), reduced chi^2 = {chi2:.2f}')
print(f'  ==> beta = {beta:.5f} +/- {beta_err:.5f}')

# --- beta(q) over every bin ---
bq, bb, bad = [], [], []
for qi in range(len(q)):
    gb = np.nanmean(G[:, :, qi], axis=0)
    eb = np.sqrt(np.nansum(E[:, :, qi]**2, axis=0)) / n_rep
    lvl, _, c2 = flat_level(tau, gb, eb)
    (bad if c2 > CHI2_FLAT_MAX else bq).append(qi)
    if c2 <= CHI2_FLAT_MAX:
        bb.append(lvl - 1.0)
print(f'  beta(q): {len(bq)} usable bins, {len(bad)} rejected as non-flat '
      f'(bins {[i+1 for i in bad]})')

# --- figure ---
apply_style()
fig, (axg, axq) = plt.subplots(1, 2, figsize=FIG_SIZE)
label_panels((axg, axq))

m = tau > 0
for i in range(n_rep):
    axg.plot(tau[m], G[i, m, Q_INDEX], '-', color='0.82', lw=LW_THIN, zorder=1)
axg.errorbar(tau[m], gbar[m], yerr=ebar[m], fmt='o', color='k', mfc='none',
             ms=MS_DENSE, mew=LW_THIN, capsize=1.5, elinewidth=LW_THIN,
             capthick=LW_THIN, zorder=2)
axg.plot(tau[m], intercept + slope * np.log10(tau[m]), 'r-', lw=LW_DATA, zorder=3)
axg.set_xscale('log')
axg.set_xlabel(r'Delay Time, $\tau$ (s)')
axg.set_ylabel('$g_2$')
axg.set_ylim(1.10, 1.16)
axg.set_title(rf'$Q = {q[Q_INDEX]:.5f}\ \AA^{{-1}}$')
add_minor_grid(axg)
axg.text(0.03, 0.06, rf'$\beta = {beta:.4f}$', transform=axg.transAxes,
         ha='left', va='bottom', color='r')

axq.plot(q[bq], bb, 'ko', mfc='none', ms=MS, mew=MEW)
axq.axhline(beta, color='r', ls='-', lw=LW_DATA)
axq.axvspan(0.00376, 0.00827, color='0.88', zorder=0)
axq.text(0.0058, 0.1465, 'range used\nfor the sample', ha='center', va='top', color='0.35')
axq.set_xlabel(r'$Q$ ($\AA^{-1}$)')
axq.set_ylabel(r'Contrast, $\beta$')
axq.set_ylim(0.10, 0.15)
add_minor_grid(axq)

fig.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig, 'FigureS4_Contrast.pdf')
plt.show()

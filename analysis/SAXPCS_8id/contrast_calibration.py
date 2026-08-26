"""Figure S4: measurement of the instrumental speckle contrast beta.

beta is a property of the beamline -- the transverse coherence of the incident
beam and the speckle-to-pixel size ratio -- not of the sample, so it is measured
once on a STATIC reference and then held fixed in every g2 fit.  A static sample
has no decorrelation, so its correlation function is flat at

    g2(tau) = 1 + beta.

Reference: nano-porous glass, 10 nm nominal pore diameter, 10 x 10 mm, 1.2 mm
thick (Doraglas S10-10-1200-50).  Dataset F0145, 6 C, 50 repeat acquisitions of
100,000 frames each, taken at 19:49 on 2022-03-06 -- about two hours before the
isothermal series of Figure 3 (21:57), in the same beamline configuration.

Panel (a): the 50 repeats and their average at q = 0.02067 A^-1, the bin of
highest scattered intensity, with the straight-line fit whose level gives beta.
Panel (b): beta from the same procedure at every q bin, showing that it is
essentially independent of q, so the value measured at the strongest bin
transfers to the low-q bins used for the sample.

Two bins are excluded from panel (b) and flagged in the print-out: bin 2
(q = 0.00489 A^-1) returns beta = 0.65 with a reduced chi^2 against a constant
of 2346, and bin 15 (q = 0.01954 A^-1) is also strongly non-flat.  Both are
detector/parasitic-scattering artefacts, not contrast.  Every one of the 50
repeats at the bin actually used is statistically flat (worst reduced chi^2 =
1.30), so no repeat is rejected.
"""

import glob
import os
import re
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)

GLASS = ('/home/8-id-i/2022-1/babnigg202203_nexus/reprocess_results/'
         'F0145_10nm_Glass_006C_att00_Rq0_*_results.hdf')
Q_INDEX = 15                 # 0-based; q = 0.02067 A^-1, the highest-intensity bin
CHI2_FLAT_MAX = 5.0          # a q bin whose averaged g2 is this far from flat is
                             # an artefact, not a contrast measurement
FIG_SIZE = (DOUBLE_COL, 2.9)


def load(files):
    G, E = [], []
    for f in files:
        with h5py.File(f, 'r') as hf:
            ft = float(np.asarray(hf['/entry/instrument/detector_1/frame_time'][()]).reshape(-1)[0])
            t = np.asarray(hf['/xpcs/multitau/delay_list'][()])
            tau = (t[:, 0] if t.ndim > 1 else t) * ft
            q = hf['/xpcs/qmap/dynamic_v_list_dim0'][()]
            G.append(hf['/xpcs/multitau/normalized_g2'][()])
            E.append(hf['/xpcs/multitau/normalized_g2_err'][()])
    return tau, np.asarray(q), np.array(G), np.array(E)


def flat_level(tau, g, e):
    """Weighted mean level of a flat g2, its error, and the reduced chi^2 of the
    flat hypothesis."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & np.isfinite(e) & (e > 0)
    w = 1.0 / e[ok]**2
    lvl = np.sum(w * g[ok]) / np.sum(w)
    return lvl, 1.0 / np.sqrt(np.sum(w)), np.sum(w * (g[ok] - lvl)**2) / (ok.sum() - 1)


files = sorted(glob.glob(GLASS))
tau, q, G, E = load(files)
runs = [int(re.search(r'Rq0_(\d+)_results', f).group(1)) for f in files]
print(f'{len(files)} repeats of the 10 nm glass standard, {len(q)} q bins')

# --- per-repeat check at the working bin: is any repeat an outlier? ---
worst = max(flat_level(tau, G[i, :, Q_INDEX], E[i, :, Q_INDEX])[2] for i in range(len(files)))
print(f'q bin {Q_INDEX+1} (q = {q[Q_INDEX]:.5f} A^-1): worst per-repeat '
      f'reduced chi^2 against a constant = {worst:.2f} -> no repeat rejected')

# --- averaged g2 at the working bin, and the straight-line fit ---
gbar = np.nanmean(G[:, :, Q_INDEX], axis=0)
ebar = np.sqrt(np.nansum(E[:, :, Q_INDEX]**2, axis=0)) / len(files)
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
    eb = np.sqrt(np.nansum(E[:, :, qi]**2, axis=0)) / len(files)
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
for i in range(len(files)):
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

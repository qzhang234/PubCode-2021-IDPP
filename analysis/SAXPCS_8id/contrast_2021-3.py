"""Figure S6: the speckle contrast of the 2021-3 beamtime, and why it is not
measured over the bins the paper fits.

(a) g2 of the static reference at the lowest fitted Q and inside the plateau.
    A static sample cannot relax, so both are flat and their level is 1 + beta.
(b) beta against Q across every dynamic bin.

THE POINT OF (b).  beta is a property of the optics and the detector and is
near enough constant in Q over this range, so a flat line is what a clean
reference returns.  This one does not give a flat line below 0.0085 A^-1: it
climbs from 0.1226 at 0.0036 A^-1 to 0.140-0.143 by 0.0086, a 12% rise.  The
reference is the NIST glassy carbon standard, whose surface scratches scatter
statically; that scattering adds to the denominator of g2 without fluctuating,
so it dilutes the speckle and depresses the APPARENT contrast, most strongly at
low Q where it is brightest.  The 2022-1 reference, a different sample, rises
by only 2% across the same interval, which places the effect on the standard
rather than on the beamline.

So beta is taken from the plateau and applied to every bin.  BETA_BAND ends at
0.0165 A^-1, which keeps out a detector module gap at 0.0175-0.0185, where beta
dips to 0.126, and the bin at 0.0214 containing the stuck pixel, where it reads
high.  Both are marked in (b).

A SINGLE SCALAR, not a value per bin: inside the plateau the bin-to-bin scatter
is instrumental noise, and letting it into beta would feed that noise straight
into f, which is the quantity the paper reports.

WHICH WAY THE OLD CHOICE ERRED.  Fitting with beta taken from the five fitted
bins, 0.1333, rather than from the plateau understates beta by 5%.  Because the
measured plateau of g2 is beta (1-f)^2, too small a beta returns too large an f
and so too small an arrested fraction: the published maxima move from 0.96 to
0.93 when the plateau value is used.  The ordering of f in Q is unchanged.

Reads data_2021-3/ only.
"""

import os
import sys

import numpy as np
import h5py
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MEW, LW_THIN, LW_DATA,     # noqa: E402
                              apply_style, add_minor_grid, label_panels, save_fig)
from nexus_read import read_g2                                           # noqa: E402

DATA = os.path.join(HERE, 'data_2021-3')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI = [0, 1, 2, 3, 4]                  # the bins the paper fits
BETA_BAND = (0.0085, 0.0165)          # A^-1, the flat part
GAP = (0.0170, 0.0190)                # detector module gap
STUCK_Q = 0.02144                     # the bin holding the stuck pixel


def flat_level(tau, g, e):
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * g[ok]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


with h5py.File(STANDARD, 'r') as hf:
    tau, g2, g2e, q = read_g2(hf)
lv = [flat_level(tau, g2[:, i], g2e[:, i]) for i in range(g2.shape[1])]
beta = np.array([x[0] - 1.0 for x in lv])
berr = np.array([x[1] for x in lv])
ok = np.isfinite(beta)
sel = ok & (q >= BETA_BAND[0]) & (q <= BETA_BAND[1])
w = 1.0 / berr[sel] ** 2
BETA = float(np.sum(w * beta[sel]) / np.sum(w))
print(f'beta = {BETA:.5f} over {sel.sum()} bins in {BETA_BAND[0]}-{BETA_BAND[1]} A^-1 '
      f'(spread {beta[sel].min():.4f}-{beta[sel].max():.4f})')
print('over the five fitted bins it would read '
      + ', '.join(f'{beta[i]:.4f}' for i in QI))

apply_style()
fig = plt.figure(figsize=(DOUBLE_COL, 2.55))
ax1 = fig.add_axes([0.075, 0.205, 0.38, 0.695])
ax2 = fig.add_axes([0.600, 0.205, 0.385, 0.695])

for qi, c, lab in ((QI[0], '#c0392b', f'{q[QI[0]]*1e3:.2f}'),
                   (int(np.argmin(abs(q - 0.0115))), '#1f77b4',
                    f'{q[int(np.argmin(abs(q - 0.0115)))]*1e3:.2f}')):
    m = (tau > 0) & (tau < 2) & np.isfinite(g2[:, qi]) & (g2e[:, qi] > 0)
    ax1.errorbar(tau[m], g2[m, qi] - 1.0, yerr=g2e[m, qi], fmt='o', color=c, ms=MS - 1,
                 mfc='none', mew=MEW, ls='none', elinewidth=LW_THIN, capsize=1.0,
                 capthick=LW_THIN, label=lab)
    ax1.axhline(beta[qi], color=c, lw=LW_DATA, zorder=3)
ax1.set_xscale('log')
ax1.set_xlabel(r'Delay Time, $\Delta t$ (s)')
ax1.set_ylabel(r'$g_2 - 1$')
ax1.set_ylim(0.09, 0.175)
leg = ax1.legend(title=r'$Q$ ($10^{-3}$ $\AA^{-1}$)', loc='upper right',
                 handlelength=0.9, handletextpad=0.35, borderpad=0.3, ncol=2,
                 columnspacing=0.7)
leg.get_title().set_fontsize(plt.rcParams['font.size'])
leg.get_frame().set_linewidth(LW_THIN)
add_minor_grid(ax1)

ax2.axvspan(BETA_BAND[0] * 1e3, BETA_BAND[1] * 1e3, color='crimson', alpha=0.10, lw=0)
ax2.axhline(BETA, color='crimson', lw=LW_DATA, zorder=2)
fitted = ok & (q <= q[QI[-1]] + 1e-6)
gapb = ok & (q >= GAP[0]) & (q <= GAP[1])
stuck = ok & (abs(q - STUCK_Q) < 5e-4)
rest = ok & ~fitted & ~gapb & ~stuck
for m, c, mk, lab in ((rest, '0.15', 'o', None),
                      (fitted, '#c0392b', 's', 'fitted bins'),
                      (gapb, '#1f77b4', '^', 'module gap'),
                      (stuck, '#2ca02c', 'v', 'stuck pixel')):
    if m.any():
        ax2.errorbar(q[m] * 1e3, beta[m], yerr=berr[m], fmt=mk, color=c, ms=MS,
                     mfc='none', mew=MEW, ls='none', elinewidth=LW_THIN, capsize=1.0,
                     capthick=LW_THIN, label=lab, zorder=3)
ax2.set_xlabel(r'$Q$ ($10^{-3}$ $\AA^{-1}$)')
ax2.set_ylabel(r'$\beta$')
ax2.set_xlim(2.5, 26)
ax2.set_ylim(0.115, 0.150)
leg2 = ax2.legend(loc='lower right', handlelength=0.9, handletextpad=0.35,
                  borderpad=0.3, labelspacing=0.2)
leg2.get_frame().set_linewidth(LW_THIN)
add_minor_grid(ax2)
label_panels([ax1, ax2])
save_fig(fig, 'FigureS6_Contrast.pdf')

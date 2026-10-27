"""The speckle contrast of the 2021-3 beamtime, and why it is not measured over
the bins the paper fits.

This prints beta; it draws no figure.  The contrast had its own SI figure up to
27 Oct 2026, and that figure was dropped: the dashed 1 + beta line already drawn
on every g2 panel shows where the contrast sits, so the figure only restated it.
How beta is measured is now given in words in SI Section 8.

beta is a property of the optics and the detector and is
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
0.0165 A^-1, the last bin before beta leaves the plateau: it dips to 0.126 at
0.0175-0.0185 and reads 0.1429 at 0.0214, against 0.1399 in the bins either
side.  The 0.0214 bin is the one holding the stuck pixel.  The cause of the dip
just below it is not established, and no cause is assumed here; both intervals
are simply outside the band beta is read from.

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

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from nexus_read import read_g2                                           # noqa: E402

DATA = os.path.join(HERE, 'data_2021-3')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI = [0, 1, 2, 3, 4]                  # the bins the paper fits
BETA_BAND = (0.0085, 0.0165)          # A^-1, the flat part


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
print('outside the band: ' + ', '.join(f'{q[i]:.4f} {beta[i]:.4f}'
                                       for i in range(len(q)) if ok[i] and not sel[i]))

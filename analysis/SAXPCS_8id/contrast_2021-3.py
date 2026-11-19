"""Print the speckle contrast beta of the 2021-3 beamtime; draws no figure.

A static sample cannot relax, so each of its correlation functions is flat at
1 + beta.  beta is read on the static standard B0005 (NIST glassy carbon) over
the flat band xpcs_fit.BETA_BAND, 0.0085-0.0165 A^-1, as in every fit.  The
script also prints the lower values of the five fitted bins and of the bins
outside the band: the standard's scratches scatter statically and depress the
apparent contrast at low Q.  See the repository README.  Reads data_2021-3/ only.
"""

import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
from nexus_read import read_g2                                           # noqa: E402
from xpcs_fit import BETA_BAND, flat_level                               # noqa: E402

DATA = os.path.join(HERE, 'data_2021-3')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
QI = [0, 1, 2, 3, 4]                  # the bins the paper fits

with h5py.File(STANDARD, 'r') as hf:
    tau, g2, g2e, q = read_g2(hf)
n_bins = g2.shape[1]
beta = np.zeros(n_bins)
berr = np.zeros(n_bins)
for i in range(n_bins):
    level, level_err = flat_level(tau, g2[:, i], g2e[:, i])
    beta[i] = level - 1.0
    berr[i] = level_err
ok = np.isfinite(beta)
sel = ok & (q >= BETA_BAND[0]) & (q <= BETA_BAND[1])
w = 1.0 / berr[sel] ** 2
BETA = float(np.sum(w * beta[sel]) / np.sum(w))
print(f'beta = {BETA:.5f} over {sel.sum()} bins in {BETA_BAND[0]}-{BETA_BAND[1]} A^-1 '
      f'(spread {beta[sel].min():.4f}-{beta[sel].max():.4f})')
fitted = []
for i in QI:
    fitted.append(f'{beta[i]:.4f}')
print('over the five fitted bins it would read ' + ', '.join(fitted))
outside = []
for i in range(len(q)):
    if ok[i] and not sel[i]:
        outside.append(f'{q[i]:.4f} {beta[i]:.4f}')
print('outside the band: ' + ', '.join(outside))

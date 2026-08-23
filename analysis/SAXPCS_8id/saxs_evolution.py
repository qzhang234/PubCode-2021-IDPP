"""SI figure: full SAXS 1d evolution for all B0147 files (D0138 subtracted).

Every B0147 averaged file is plotted, coloured by elapsed time (time origin =
first B0147 file, frames 1-200).  B0146 (6 C reference, before the 30 C
isothermal) is shown as a grey dashed curve.  The legend is placed outside, on
the right-hand side of the axes.
"""

import glob
import os
import re
import sys
from datetime import datetime

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Absolute scattering cross-section calibration lives in the shared module so
# this figure uses exactly the same constants and per-file coefficient as the
# main saxpcs.py figure.  See abs_xsec.py for the full derivation.
from abs_xsec import abs_xsec_coef, calibration_summary

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.prl_style import SINGLE_COL, apply_style, add_minor_grid, save_tight

# --- PARAMETERS ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")
BACKGROUND_HEADER = 'D0138'
# Absolute scattering cross-section: convert the raw SAXS 1d to an absolute
# differential cross section (cm^-1) as
#     I_abs(q) = coef_sam * I_sample(q) - coef_buf * I_buffer(q)
# coef_sam / coef_buf are NOT hard-coded any more: each is computed by
# abs_xsec_coef() from that file's own range-averaged ion-chamber readings, so a
# drift in incident flux across the time series is handled per file.  coef_buf
# is the D0138 buffer's own coefficient; coef_sam is recomputed for every curve.
PHI_AVERAGE = True
CMAP = plt.cm.turbo

SAXS_PATH       = '/xpcs/temporal_mean/scattering_1d'
STATIC_MAP_PATH = '/xpcs/qmap/static_index_mapping'
STATIC_Q_PATH   = '/xpcs/qmap/static_v_list_dim0'
STATIC_PHI_PATH = '/xpcs/qmap/static_v_list_dim1'
START_TIME_PATH = '/entry/start_time'
TIME_FORMAT     = '%Y-%m-%d %H:%M:%S'

_name_re = re.compile(r'Average_([A-Za-z]\d+)_.*?_(\d+)_(\d+)_results')


def parse_name(fname):
    m = _name_re.search(os.path.basename(fname))
    return (m.group(1), int(m.group(2)), int(m.group(3))) if m else (None, -1, -1)


def read_start_time(hf):
    raw = hf[START_TIME_PATH][()]
    if isinstance(raw, np.ndarray):
        raw = raw.reshape(-1)[0]
    if isinstance(raw, bytes):
        raw = raw.decode('utf-8')
    return datetime.strptime(str(raw).strip(), TIME_FORMAT)


def read_saxs_iq(hf, phi_average=True):
    intensity = np.asarray(hf[SAXS_PATH][()]).reshape(-1)
    idx_map = hf[STATIC_MAP_PATH][()]
    q_list = hf[STATIC_Q_PATH][()]
    n_phi = hf[STATIC_PHI_PATH].shape[0]
    q_idx = idx_map // n_phi
    uq = np.unique(q_idx)
    if phi_average and n_phi > 1:
        inten = np.array([np.nanmean(intensity[q_idx == qi]) for qi in uq])
    else:
        inten = np.array([intensity[q_idx == qi][0] for qi in uq])
    return q_list[uq], inten


# --- DISCOVER FILES ---
file_paths = sorted(glob.glob(os.path.join(data_dir, '*.hdf')))
by_header, start_times = {}, {}
for fp in file_paths:
    header = parse_name(fp)[0]
    with h5py.File(fp, 'r') as hf:
        start_times[fp] = read_start_time(hf)
    by_header.setdefault(header, []).append(fp)
for h in by_header:
    by_header[h].sort(key=lambda p: parse_name(p)[1])

b0147 = by_header['B0147']
t_ref = start_times[b0147[0]]
elapsed = {fp: (start_times[fp] - t_ref).total_seconds() for fp in b0147}
norm = Normalize(vmin=0, vmax=max(elapsed.values()))

print(calibration_summary())

# Buffer curve AND its own absolute coefficient (coef_buf), both from D0138.
bg_I = None
coef_buf = None
if BACKGROUND_HEADER in by_header:
    with h5py.File(by_header[BACKGROUND_HEADER][0], 'r') as hf:
        _, bg_I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_buf = abs_xsec_coef(hf)
    print(f'coef_buf ({BACKGROUND_HEADER}) = {coef_buf:.3e}')

# --- FIGURE ---
apply_style()
fig, ax = plt.subplots(figsize=(SINGLE_COL, 3.3))

for fp in by_header.get('B0146', []):
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                       # this file's own coefficient
    I = coef_sam * I - coef_buf * bg_I if bg_I is not None else coef_sam * I
    pos = I > 0
    ax.plot(q[pos], I[pos], color='0.5', marker='s', ls='none', ms=3,
            mfc='none', mew=0.5, label=r'6$^{\circ}$C Ref')

for fp in b0147:
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                       # this file's own coefficient
    I = coef_sam * I - coef_buf * bg_I if bg_I is not None else coef_sam * I
    pos = I > 0
    ax.plot(q[pos], I[pos], color=CMAP(norm(elapsed[fp])), marker='o', ls='none',
            ms=3, mfc='none', mew=0.5, label=f'{elapsed[fp]:.0f} s')

ax.set_xscale('log')
ax.set_yscale('log')
ax.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax.set_ylabel(r'$d\Sigma/d\Omega$ (cm$^{-1}$)')
add_minor_grid(ax)
ax.legend(loc='upper right', title='Elapsed time',
          ncol=2, columnspacing=1.0, handletextpad=0.4, labelspacing=0.25)

save_tight(fig, 'FigureS6_SAXS_Evolution.pdf')
plt.show()
print('wrote FigureS6_SAXS_Evolution.pdf')

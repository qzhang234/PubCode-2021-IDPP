"""The waiting-time clock and the absolute SAXS of the 2021-3 temperature series.

Shared by temperature_series_2021-3.py, arrest_vs_temperature.py and
general_model_check.py, so all of them use one t_w and one I(Q):

- t_w = 0 is the first acquisition with the stage RTD within AT_TARGET of its
  target, and each averaged group sits at the midpoint of its acquisitions;
- I(Q) is absolute (abs_xsec.py) with BG_COEF times the empty cell E0411
  subtracted.

The repository README explains both choices.  Reads data_2021-3/ only.
"""

import csv
import glob
import os
from datetime import datetime

import h5py
import numpy as np

from abs_xsec import abs_xsec_coef
from nexus_read import parse_name, read_iq

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data_2021-3')
TRACE = os.path.join(DATA, 'trace_2021-3.csv')
BACKGROUND = os.path.join(DATA, 'Average_E0411_Empty_3_006C_att00_Rq0_00001_00050_results.hdf')
BAND = (0.004, 0.008)                 # A^-1, the low-Q band of Figure 3c and of the trace
AT_TARGET = 0.01                      # C
BG_COEF = 0.80                        # fraction of the empty cell subtracted (README)
# One hue and one marker shape per temperature, 28.00 to 28.75 C: blue, green,
# orange, red (through acs_style.hsv_colors); circle, square, triangle, diamond.
T_HUES = (0.66, 0.33, 0.09, 0.0)
T_MARKERS = ('o', 's', '^', 'D')


# --- the clock ---------------------------------------------------------------
def trace_rows(run):
    """The rows of trace_2021-3.csv for one run, in acquisition order."""
    rows = []
    with open(TRACE) as fh:
        for row in csv.DictReader(fh):
            if row['run'] == run:
                rows.append(row)
    frames = [int(row['frame']) for row in rows]
    order = np.argsort(frames, kind='stable')
    return [rows[k] for k in order]


def group_files(run):
    """The averaged files of one run, in order of their first acquisition."""
    files = glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf'))
    first_frames = [parse_name(fp)[1] for fp in files]
    order = np.argsort(first_frames, kind='stable')
    return [files[k] for k in order]


def ramp_reference(run):
    """The ~20 C reference of one run: its first 20 acquisitions (RTD 20.0-20.2 C),
    averaged by reduce_babnigg202110.py ramp20 into data_2021-3/ramp20/."""
    files = glob.glob(os.path.join(DATA, 'ramp20', f'Average_{run}_*_results.hdf'))
    if len(files) != 1:
        raise FileNotFoundError(f'{len(files)} ~20 C references for {run} in data_2021-3/ramp20/')
    return files[0]


def acquisitions(run):
    """Frame numbers, seconds from t_w = 0, RTD readings, and the t_w = 0 datetime."""
    rows = trace_rows(run)
    frame = np.array([int(r['frame']) for r in rows])
    when = [datetime.strptime(r['start_time'], '%Y-%m-%d %H:%M:%S') for r in rows]
    T = np.array([float(r['qnw3_C']) for r in rows])
    on_target = T >= float(rows[0]['qnw3_set_C']) - AT_TARGET
    if not on_target.any():
        raise ValueError(f'{run} never reaches its target temperature, so t_w = 0 is undefined')
    i0 = int(np.argmax(on_target))
    t = np.array([(w - when[i0]).total_seconds() for w in when])
    return frame, t, T, when[i0]


def group_midpoint(fp, frame, t):
    """Midpoint, in s from t_w = 0, of the acquisitions an averaged file contains."""
    header, first, last = parse_name(fp)
    inside = (frame >= first) & (frame <= last)
    return float(0.5 * (t[inside][0] + t[inside][-1]))


# --- the SAXS ----------------------------------------------------------------
# the empty cell, absolute
with h5py.File(BACKGROUND, 'r') as hf:
    Q_BG, I_BG = read_iq(hf)
    COEF_BG = abs_xsec_coef(hf)


def abs_iq(fp):
    """Absolute, background-subtracted I(Q) of one averaged file, in mm^-1."""
    with h5py.File(fp, 'r') as hf:
        q, I = read_iq(hf)
        coef = abs_xsec_coef(hf)
    if not np.allclose(q, Q_BG):
        raise ValueError(f'{fp} and the background are on different Q grids')
    return q, coef * I - BG_COEF * COEF_BG * I_BG


def band_intensity(fp):
    """Mean absolute I(Q) over BAND of one averaged file, in mm^-1."""
    q, I = abs_iq(fp)
    return float(np.mean(I[(q >= BAND[0]) & (q < BAND[1])]))

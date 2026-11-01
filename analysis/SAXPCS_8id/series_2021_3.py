"""Shared pieces of the 2021-3 temperature series: the clock and the SAXS.

temperature_series_2021-3.py (Figures 2, S4-S7 and the fast-mode SI figure)
and arrest_vs_temperature.py (Figure 3 and the SI collapse figure) need the
same waiting time and the same absolute I(Q), so both are defined here once.

THE CLOCK.  t_w = 0 is the first acquisition with the stage RTD at its target.
A group averages about 100 consecutive acquisitions spread over about 10 min,
so f and I(Q) of a group are averages over that window, and the group is
placed at the MIDPOINT of the acquisitions it contains.  Placing it at its
first acquisition instead puts the first group of a hold at 6 s rather than
near 300 s, a decade early on the logarithmic axes of the scaling analysis,
and it biases T_s low, because a fixed offset of about 300 s is a larger
fraction of the short arrest times of the warmer holds.

THE SAXS.  Absolute, background-subtracted I(Q): the empty cell E0411
subtracted with the coefficient BG_COEF, both explained below.

Reads data_2021-3/ only.
"""

import csv
import os
from datetime import datetime

import h5py
import numpy as np

from abs_xsec import abs_xsec_coef
from nexus_read import parse_name

HERE = os.path.dirname(os.path.abspath(__file__))
DATA = os.path.join(HERE, 'data_2021-3')
TRACE = os.path.join(DATA, 'trace_2021-3.csv')
BACKGROUND = os.path.join(DATA, 'Average_E0411_Empty_3_006C_att00_Rq0_00001_00050_results.hdf')
BAND = (0.004, 0.008)                 # A^-1, the low-Q band of the trace and of Figure 3c
AT_TARGET = 0.01                      # C
# One hue and one marker shape per hold, 28.00 to 28.75 C, cold to warm: blue,
# green, orange, red; circle, square, triangle, diamond.  Both figure scripts
# mark the holds with these (the hues through acs_style.hsv_colors), so a
# temperature looks the same in every figure and never rests on colour alone.
T_HUES = (0.66, 0.33, 0.09, 0.0)
T_MARKERS = ('o', 's', '^', 'D')


# --- the clock ---------------------------------------------------------------
def acquisitions(run):
    """Frame numbers, seconds from t_w = 0, RTD readings, and the t_w = 0 datetime."""
    rows = sorted((r for r in csv.DictReader(open(TRACE)) if r['run'] == run),
                  key=lambda r: int(r['frame']))
    frame = np.array([int(r['frame']) for r in rows])
    when = [datetime.strptime(r['start_time'], '%Y-%m-%d %H:%M:%S') for r in rows]
    T = np.array([float(r['qnw3_C']) for r in rows])
    i0 = int(np.argmax(T >= float(rows[0]['qnw3_set_C']) - AT_TARGET))
    t = np.array([(w - when[i0]).total_seconds() for w in when])
    return frame, t, T, when[i0]


def group_midpoint(fp, frame, t):
    """Midpoint, in s from t_w = 0, of the acquisitions an averaged file contains."""
    _, a, b = parse_name(fp)
    m = (frame >= a) & (frame <= b)
    return float(0.5 * (t[m][0] + t[m][-1]))


# --- the SAXS ----------------------------------------------------------------
def read_iq(hf):
    """Azimuthally averaged I(Q), as nexus_read.read_saxs_iq computes it.

    One value per (Q, phi) bin, flattened, phi running fastest, averaged over
    phi at each Q.
    """
    intensity = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()], dtype=float).reshape(-1)
    idx = hf['/xpcs/qmap/static_index_mapping'][()]
    q_list = hf['/xpcs/qmap/static_v_list_dim0'][()]
    n_phi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
    q_idx = idx // n_phi
    uq = np.unique(q_idx)
    return q_list[uq], np.array([np.nanmean(intensity[q_idx == i]) for i in uq])


with h5py.File(BACKGROUND, 'r') as _hf:
    Q_BG, I_BG = read_iq(_hf)
    COEF_BG = abs_xsec_coef(_hf)

# --- THE CHOSEN COEFFICIENT, AND WHY IT IS NOT THE FLATTEST ONE -------------
# BG_COEF slides along a single trade-off that this dataset cannot escape, and
# 0.80 is the agreed point on it.  Measured on the four 6 C references:
#
#     BG_COEF   rise (low Q / high Q)   worst upward step
#      0.90           0.86-0.98                9-14%
#      0.865          0.94-1.00                6-8%
#      0.80           1.24-1.44                1-4%     <-- used
#      0.70           1.62-1.87                0-1%
#      0.60           1.98-2.27                0%
#
# The flattest coefficient is 0.865, but it leaves a 6-8% dip near
# 0.004-0.005 A^-1, and a dip is a subtraction artefact rather than anything
# the sample does.  Driving the dip out entirely takes 0.60, which buys a
# factor of two climb toward low Q instead.  0.80 keeps the climb under 1.5
# and the dip inside a few per cent, which is the scatter of the points
# themselves.
#
# WHAT NO COEFFICIENT BUYS.  The 2022-1 buffer-subtracted 6 C reference fell
# 3.25-fold from 0.004 to 0.029 A^-1 with a log-log curvature of -0.46: slopes
# STEEPEN from -0.49 to -1.07, the roll-over of a form factor.  The 2021-3
# reference has the opposite curvature at every coefficient, a steep low-Q tail
# standing on a plateau.  A BUFFER measured in the same cell takes windows,
# solvent and parasitic scattering off together; 2021-3 has no buffer, and at
# 0.004 A^-1 its empty cell is 3.0 to 3.7 times the 6 C signal (2.4 to 3.0
# after BG_COEF), so subtracting enough to remove the instrumental tail removes
# the sample's own low-Q signal with it.  That is a limit of the 2021-3
# measurement, not of the coefficient, and it is confined to the 6 C
# reference: every hold group averages 54 to 322 times the background over
# 0.004-0.008 A^-1.  temperature_series_2021-3.py prints the flattest
# coefficient of each reference, the rise and step at BG_COEF, and both ratios.
BG_COEF = 0.80


def abs_iq(fp):
    """Absolute, background-subtracted I(Q) of one averaged file, in mm^-1."""
    with h5py.File(fp, 'r') as hf:
        q, I = read_iq(hf)
        coef = abs_xsec_coef(hf)
    assert np.allclose(q, Q_BG), 'sample and background on different Q grids'
    return q, coef * I - BG_COEF * COEF_BG * I_BG


def band_intensity(fp):
    """Mean absolute I(Q) over BAND of one averaged file, in mm^-1."""
    q, I = abs_iq(fp)
    return float(np.mean(I[(q >= BAND[0]) & (q < BAND[1])]))

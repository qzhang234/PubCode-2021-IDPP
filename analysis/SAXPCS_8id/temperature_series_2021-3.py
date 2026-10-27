"""Figures S11-S14: the 2021-3 temperature series, drawn as Figure 3 and Figure S9 are.

Reads only data_2021-3/, written by reduce_babnigg202110.py.  One figure per
temperature, in the order 28.00, 28.25, 28.50, 28.75 C (loadings D8, D6, D7,
D9).  Each figure carries the three panels of Figure 3 on its top row and the
four further Q bins of Figure S9 on its bottom row:

  (a) absolute, background-subtracted I(Q): the 6 C reference of that loading,
      the first group on target (t_w = 0) and the groups shown in (b);
  (b) g2 at the lowest Q with the fits;
  (c) fast fraction f against t_w in the five Q bins, with the fitted
      diffusion coefficient D0 of the fast mode on the right axis;
  (d)-(g) g2 and fits at the four other Q bins.

THE DATA.  Four loadings in one cell block, each taken from 20 C to its final
temperature at 0.10 C/min and then held there; temperatures are the RTD
readings of the QNW stage.  t_w = 0 is the first acquisition at the target,
the analogue of the RTD reaching 30 C in Figure 3.

WHAT IS THE SAME AS FIGURE 3.  The reduction (groups of 100 acquisitions through
the two outlier cuts of average_ranges.py), the model (Equation 2, the fast mode
tied to 1/(D0 Q^2), beta fixed), the five Q bins, the error-weighted least
squares and covariance, the chi^2 guard of 10, and the absolute scale: every
I(Q) is converted with abs_xsec.py from its own ion-chamber readings, which
2021-3 recorded on the same chambers, and the background is subtracted as
coef_sample * I_sample - coef_bg * I_bg, exactly as Figure 3 subtracts the buffer.

WHAT DIFFERS, AND WHY.
* The background is E0411.  Six candidates exist; E0411, G0409, G0410 and
  E0408 were recorded on 17 Oct, the same day as the samples.  THE NAMES ARE
  NOT EVIDENCE: they were entered at the beamline and do not reliably
  separate a cell left empty from one holding buffer, so the choice rests on
  the monitors and on I(Q).  Four of the six recorded a transmitted beam.
  E0411 passes 0.94 and the other three pass 1.00, the open-path value, so
  only E0411 has material in the beam for the chambers to see.  None passes
  anything near 0.41, which is what a 2 mm buffer column behind the same
  material would give, so none of them holds a liquid column of the sample
  thickness.  Fitting the three that pass 1.00 to a 6 C reference needs a
  coefficient of 3 to 5 and still leaves 7-28% against 3-5% for E0411, and
  averaging them does not help, so what they lack is shape, not statistics.
  Had any of them held buffer it would matter little over this Q range: the
  absorption goes into the coefficient and 2 mm of water scatters about
  0.0016 mm^-1, 6% of the plateau and flat.
* The empty cell is subtracted with a coefficient BG_COEF = 0.80.  The
  flattest coefficient is 0.865, but it leaves a dip near 0.004-0.005 A^-1,
  and a dip is a subtraction artefact rather than anything the sample does;
  0.80 is the agreed point on that trade-off.  Subtracting all of the empty
  cell leaves a bowl of 40-60% near 0.004 A^-1, deepest where the background
  is steepest.  The comment at the fit tabulates the trade-off, gives the
  estimator and says why ordinary least squares on I gets this wrong.
* Panel (a) starts at the lowest fitted Q bin.  Below it the points sit on the
  edge of the beamstop shadow, where no single coefficient reconciles the empty
  cell with the sample and the subtraction leaves a spike.
* D8 has two 6 C measurements and Figure S11 shows the LATER one, H0403.
  H0402, taken before the hold, was the first measurement of the beamtime and
  is not a dissolved state: it carries a Q^-2.7 low-Q excess that dies out
  above 0.01 A^-1, where it is 4% of the plateau, and it fits a background
  coefficient of 1.15 against 0.85-0.88 for every other reference.  H0403 is
  the same loading in the same cell, returned to 6 C after the 28.00 C hold.
  It is flat and fits 0.860, so the excess is undissolved material in the
  earlier measurement rather than anything in the subtraction, and H0403 is
  the measurement that shows what D8 looks like dissolved.
* Pixel (y, x) = (63, 166) was stuck at 1.5 photons per frame in the 17 Oct
  measurements.  It is masked in the qmap this series is correlated against,
  ..._Lin_CrossRm2.hdf, so nothing downstream has to work around it.  It sat in
  the static bin at Q = 0.0214 A^-1, phi sector 7, and in dynamic bin 18, far
  from the five fitted bins, and masking it left every fitted number unmoved.
* beta is measured for this beamtime on the static standard B0005, but over
  0.0085-0.0165 A^-1 rather than over the five fitted bins: that standard is
  scratched and its parasitic scattering depresses the apparent contrast by 12%
  at the lowest fitted Q.  See the comment at BETA_BAND.
* p_fast is fitted in each group, started from the minimum of a p_fast profile.
  The holds run from a mobile state (p_fast near 0.9) into arrest (p_fast near
  0.3), and a single shared value is rejected by an F test in every loading.
* The rate is reported as D0 itself, the fitted coefficient of
  tau_fast = 1/(D0 Q^2).  It is a parameter of the fit and needs no further
  construction; no mean or median of the stretched decay is formed, because
  either one is dominated by the tail and swings with p_fast.
* A group is fitted when its correlation function is no noisier than the
  noisiest group Figure 3 fits: median g2_err at the lowest Q over 1e-4 to
  1e-2 s of at most 0.016.  Where more than seven groups qualify, (a), (b) and
  (d)-(g) show every second one; (c) shows them all.
"""

import csv
import glob
import os
import sys
from datetime import datetime

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from matplotlib.ticker import LogLocator, FuncFormatter

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
sys.path.insert(0, os.path.join(HERE, '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels,
                              q_log_ticks, save_fig)
import xpcs_fit                                                  # noqa: E402
from nexus_read import parse_name, read_start_time, read_g2      # noqa: E402
from abs_xsec import (abs_xsec_coef, CAL_A, CAL_B, T_EXP,         # noqa: E402
                      AIR_TRANSMISSION, NUM_FRAMES, SAMPLE_THICKNESS,
                      DELTA_OMEGA)

DATA = os.path.join(HERE, 'data_2021-3')
TRACE = os.path.join(DATA, 'trace_2021-3.csv')
STANDARD = os.path.join(DATA, 'Average_B0005_12IDB_Standard_006C_att00_Rq0_00001_00100_results.hdf')
BACKGROUND = os.path.join(DATA, 'Average_E0411_Empty_3_006C_att00_Rq0_00001_00050_results.hdf')
BAND = (0.004, 0.008)                 # A^-1, the band of Figure S4b and of the trace

# header -> (run name in the trace, sample, target C), in order of temperature
RUNS = {'H0402': ('H0402_D8_280C10p', 'D8', 28.00),
        'H0407': ('H0407_D6_282C10p', 'D6', 28.25),
        'H0404': ('H0404_D7_285C10p', 'D7', 28.50),
        'H0413': ('H0413_D9_287C10p', 'D9', 28.75)}
REF_6C = os.path.join(DATA, 'Average_H0404_D7_060C10p_att00_Rq0_00001_00020_results.hdf')
QI = [0, 1, 2, 3, 4]                  # the five Q bins, as in Figure 3
Q_MARKERS = ['o', 's', '^', 'D', 'v']
T_MARKERS = ['o', 's', '^', 'D']
T_CMAP, T_LO, T_HI = plt.cm.viridis, 0.05, 0.85
AT_TARGET = 0.01                      # C
ERR_MAX = 0.016                       # the noisiest group Figure 3 fits
CHI2_MAX = 10.0                       # the Figure 3 guard
P_GRID = np.round(np.arange(0.20, 1.31, 0.05), 2)
SHOW_RUN = 'H0404'                    # the run shown in panels (a) and (d)
COLOR_6C = '#1f77b4'                  # the 6 C reference, as in Figure 3a
FIG_SIZE = (DOUBLE_COL, 5.0)

T_COLOR = {h: T_CMAP(c) for h, c in zip(RUNS, np.linspace(T_LO, T_HI, len(RUNS)))}
INC = '/entry/instrument/incident_beam/incident_beam_intensity'
TRN = '/entry/instrument/incident_beam/transmitted_beam_intensity'


# --- beta for this beamtime -------------------------------------------------
# NOT measured over the five fitted bins.  The 2021-3 static reference is the
# NIST glassy carbon standard, whose surface scratches add parasitic scattering
# at low Q.  That scattering is static, so it dilutes the speckle and depresses
# the APPARENT contrast exactly where the paper fits.  The measurement shows it
# plainly: beta climbs from 0.1226 at 0.0036 A^-1 to a plateau of 0.140-0.143
# by 0.0086 A^-1, a 12% rise.  The 2022-1 reference, a different sample, rises
# by only 2% over the same interval, so the effect belongs to the standard and
# not to the beamline.
#
# Contrast is a property of the optics and the detector, near enough constant
# in Q over this range, so it is taken from the clean plateau and applied to
# every bin.  BETA_BAND stops at 0.0165 A^-1, which excludes a detector module
# gap at 0.0175-0.0185 where beta dips to 0.126.  A single scalar rather than a
# per-bin value: the
# bin-to-bin scatter inside the plateau is instrumental noise, and letting it
# into f would inject that noise into the very ordering the paper reports.
BETA_BAND = (0.0085, 0.0165)          # A^-1, the flat part of the standard


def flat_level(tau, g, e):
    """Inverse-variance weighted mean of g2 over the measurable delays, and its
    standard error.  Every delay carries its own g2 uncertainty, so each enters
    with weight 1/e^2 and the noisy short delays count for little."""
    ok = (tau > 0) & (tau < 2) & np.isfinite(g) & (e > 0)
    w = 1.0 / e[ok] ** 2
    return float(np.sum(w * g[ok]) / np.sum(w)), float(np.sqrt(1.0 / np.sum(w)))


def measure_beta(path, label):
    with h5py.File(path, 'r') as hf:
        tau, g2, g2e, q_std = read_g2(hf)
    lv = [flat_level(tau, g2[:, i], g2e[:, i]) for i in range(g2.shape[1])]
    b = np.array([x[0] - 1.0 for x in lv])
    e = np.array([x[1] for x in lv])
    sel = (q_std >= BETA_BAND[0]) & (q_std <= BETA_BAND[1]) & np.isfinite(b)
    w = 1.0 / e[sel] ** 2
    beta = float(np.sum(w * b[sel]) / np.sum(w))
    print(f'beta ({label}) = {beta:.5f} from {sel.sum()} bins over '
          f'{BETA_BAND[0]}-{BETA_BAND[1]} A^-1 (spread {b[sel].min():.4f}-{b[sel].max():.4f})')
    print('   over the five FITTED bins it would read '
          + ', '.join(f'{b[i]:.4f}' for i in QI) + '  <- depressed by the standard')
    return beta


BETA = measure_beta(STANDARD, '2021-3 static standard')
xpcs_fit.CONTRAST = BETA              # double_exp reads it at call time


# --- absolute scale ---
def coef_from(up_ic, dn_ic):
    """abs_xsec_coef for monitor readings given as numbers rather than an open file."""
    flux = CAL_A * (up_ic / T_EXP) + CAL_B
    transmission = dn_ic / up_ic / AIR_TRANSMISSION
    return NUM_FRAMES / (T_EXP * flux * transmission * SAMPLE_THICKNESS * DELTA_OMEGA)


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


def band_mean(hf):
    """Mean of the stored (Q, phi) entries in BAND, as the trace records it."""
    s = np.asarray(hf['/xpcs/temporal_mean/scattering_1d'][()]).reshape(-1)
    idx = hf['/xpcs/qmap/static_index_mapping'][()]
    q = hf['/xpcs/qmap/static_v_list_dim0'][()][idx // hf['/xpcs/qmap/static_v_list_dim1'].shape[0]]
    return float(np.nanmean(s[(q >= BAND[0]) & (q < BAND[1])]))


with h5py.File(BACKGROUND, 'r') as hf:
    Q_BG, I_BG = read_iq(hf)
    COEF_BG = abs_xsec_coef(hf)
    BAND_BG = band_mean(hf)
    check = coef_from(float(np.asarray(hf[INC][()]).ravel()[0]),
                      float(np.asarray(hf[TRN][()]).ravel()[0]))
assert abs(check / COEF_BG - 1) < 1e-12, 'coef_from disagrees with abs_xsec_coef'
print(f'background E0411: coefficient {COEF_BG:.4g}, absolute band mean '
      f'{COEF_BG * BAND_BG:.4f} mm^-1')


# --- how much of the background to subtract ---
# At 6 C the sample is dissolved, so its background-subtracted I(Q) has to be
# flat across the plotted range.  That fixes one number, ALPHA, in
#     I_6C(Q) - ALPHA * I_bg(Q) = plateau
#
# HOW IT IS FITTED, AND WHY NOT BY ORDINARY LEAST SQUARES.  Straight least
# squares on I minimises ABSOLUTE residuals, and at the bottom of the range the
# background is six times the 6 C signal, so the two or three lowest Q points
# carry almost all the weight and drag ALPHA up.  That is what left a 10-20%
# bowl near 0.004-0.005 A^-1 in every loading: too much background subtracted
# over the part of the range the fit was not really looking at.  Instead the
# profile is reduced to the median in each of FIT_NBINS logarithmic Q bins,
# divided by its own plateau, and ALPHA is the value on ALPHA_GRID that
# minimises the rms departure of those bin values from 1.  Every decade of Q
# then counts equally and the bins average down the point scatter.  The
# residual falls from 12-20% to 3-5% and ALPHA from about 0.89 to 0.865.
#
# ALPHA is the MEDIAN over the four loadings, one shared number for all four.
# Giving each loading its own value changes the flatness by well under a per
# cent, so the shared value costs nothing and leaves one number to report.
#
# What ALPHA cannot fix: below about 0.005 A^-1 the background is several times
# the 6 C signal, so a few per cent of shape mismatch between it
# and the sample cell is tens of per cent of the dissolved-state intensity, and
# no single coefficient removes it.  A q-scale mismatch between the two was
# tested and rejected: flattening the residual that way needs a 6% stretch of
# the background q axis and improves the rms by half a point.  The arrested and
# assembling groups are 40 to 800 times the background over the same range, so
# this limit touches the 6 C reference alone.
FIT_BAND = (0.0036, 0.0305)           # A^-1; starts at the lowest fitted Q bin
FIT_NBINS = 9                         # logarithmic Q bins across FIT_BAND
PLATEAU_BAND = (0.012, 0.0305)        # A^-1, where every reference is flat
ALPHA_GRID = np.arange(0.70, 1.15, 0.0025)


def flatness(q, y):
    """RMS fractional departure of the binned median from the plateau."""
    plateau = np.median(y[(q >= PLATEAU_BAND[0]) & (q < PLATEAU_BAND[1])])
    edges = np.exp(np.linspace(np.log(FIT_BAND[0]), np.log(FIT_BAND[1]), FIT_NBINS + 1))
    binned = np.array([np.median(y[(q >= a) & (q < b)]) for a, b in zip(edges[:-1], edges[1:])])
    return float(np.sqrt(np.mean((binned / plateau - 1.0) ** 2)))


def shape(q, y):
    """(rise, worst upward step) of the binned median across FIT_BAND.

    rise  = lowest bin / highest bin.  1 is flat; above 1 the profile climbs
            toward low Q.
    step  = largest fractional INCREASE from one bin to the next going up in Q.
            0 means the profile never turns back up, i.e. no oscillation.
    """
    edges = np.exp(np.linspace(np.log(FIT_BAND[0]), np.log(FIT_BAND[1]), FIT_NBINS + 1))
    binned = np.array([np.median(y[(q >= a) & (q < b)]) for a, b in zip(edges[:-1], edges[1:])])
    step = float(max(0.0, (np.diff(binned) / binned[:-1]).max()))
    return float(binned[0] / binned[-1]), step


def flattest_coef(ref_files):
    """Per-loading coefficient that makes the 6 C reference flattest."""
    out = {}
    for name, fp in ref_files.items():
        with h5py.File(fp, 'r') as hf:
            I = abs_xsec_coef(hf) * read_iq(hf)[1]
        cost = [flatness(Q_BG, I - a * COEF_BG * I_BG) for a in ALPHA_GRID]
        k = int(np.argmin(cost))
        out[name] = (float(ALPHA_GRID[k]), cost[k])
    return out


# The 6 C reference of each loading.  D8's is H0403, taken at 6 C after its
# 28.00 C hold, rather than H0402 taken before it: H0402 was the first
# measurement of the beamtime and still carries undissolved material, a Q^-2.7
# low-Q excess that dies out above 0.01 A^-1.  It is the same loading in the
# same cell, and it is the only one of the five 6 C measurements that is not a
# dissolved state, so H0403 is the one that shows what D8 looks like dissolved.
REFS = {'H0402': 'H0403_D8_060C10p', 'H0407': 'H0407_D6_060C10p',
        'H0404': 'H0404_D7_060C10p', 'H0413': 'H0413_D9_060C10p'}

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
# WHAT NO COEFFICIENT BUYS.  Figure 3a's 6 C reference falls 3.25-fold from
# 0.004 to 0.029 A^-1 with a log-log curvature of -0.46: slopes STEEPEN from
# -0.49 to -1.07, the roll-over of a form factor.  The 2021-3 reference has
# the opposite curvature at every coefficient, a steep low-Q tail standing on
# a plateau.  Figure 3a subtracts a BUFFER measured in the same cell, which
# takes windows, solvent and parasitic scattering off together.  2021-3 has
# no buffer, and at 0.004 A^-1 its empty cell is 4.7 times the 6 C signal, so
# subtracting enough to remove the instrumental tail removes the sample's own
# low-Q signal with it.  That is a limit of the 2021-3 measurement, not of
# the coefficient, and it is confined to the 6 C reference: every hold group
# is 40 to 800 times the background over the same range.
BG_COEF = 0.80

_BG_FITS = flattest_coef({h: os.path.join(DATA, f'Average_{r}_att00_Rq0_00001_00020_results.hdf')
                          for h, r in REFS.items()})
print(f'background coefficient: using {BG_COEF:.3f}; flattest would be '
      + ', '.join(f'{h} {a:.3f}' for h, (a, _) in _BG_FITS.items()))
for h, r in REFS.items():
    with h5py.File(os.path.join(DATA, f'Average_{r}_att00_Rq0_00001_00020_results.hdf'), 'r') as hf:
        I = abs_xsec_coef(hf) * read_iq(hf)[1]
    rise, step = shape(Q_BG, I - BG_COEF * COEF_BG * I_BG)
    print(f'   {r}: at {BG_COEF:.2f} the 6 C reference rises {rise:.2f}x toward low Q '
          f'and turns back up by at most {step * 100:.0f}%')


def abs_iq(fp):
    """Absolute, background-subtracted I(Q) of one averaged file, in mm^-1."""
    with h5py.File(fp, 'r') as hf:
        q, I = read_iq(hf)
        coef = abs_xsec_coef(hf)
    assert np.allclose(q, Q_BG), 'sample and background on different Q grids'
    return q, coef * I - BG_COEF * COEF_BG * I_BG


# --- per-acquisition trace: absolute band intensity, temperature, time ---
def load_trace(run):
    rows = [r for r in csv.DictReader(open(TRACE)) if r['run'] == run]
    rows.sort(key=lambda r: int(r['frame']))
    when = [datetime.strptime(r['start_time'], '%Y-%m-%d %H:%M:%S') for r in rows]
    T = np.array([float(r['qnw3_C']) for r in rows])
    target = float(rows[0]['qnw3_set_C'])
    I = np.array([coef_from(float(r['incident']), float(r['transmitted'])) * float(r['I_band'])
                  for r in rows]) - BG_COEF * COEF_BG * BAND_BG
    i0 = int(np.argmax(T >= target - AT_TARGET))
    t = np.array([(w - when[i0]).total_seconds() for w in when])
    return t, T, I, when[i0]


def binned(x, y, edges):
    """Median of y in each bin of x, for bins holding at least five points."""
    xc, yc = [], []
    for lo, hi in zip(edges[:-1], edges[1:]):
        m = (x >= lo) & (x < hi)
        if m.sum() >= 5:
            xc.append(np.median(x[m]))
            yc.append(np.median(y[m]))
    return np.array(xc), np.array(yc)


# --- Equation 2 in one block, p_fast free, started at its profile minimum ---
def fit_block(tau, g2, g2e, q):
    prof = [(r['red_chi2'], pf, r) for pf in P_GRID
            for r in [xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q, p1_fixed=pf)]
            if r is not None]
    chi0, pf0, best = min(prof, key=lambda x: x[0])
    xpcs_fit.P_EXP_P0[0] = pf0
    free = xpcs_fit.fit_g2_diffusive(tau, g2, g2e, QI, q)
    xpcs_fit.P_EXP_P0[0] = 0.5
    if free is not None and free['red_chi2'] <= chi0:
        best = free
    return best


traces, series = {}, {}
for h, (run, sample, target) in RUNS.items():
    t, T, I, t0 = load_trace(run)
    traces[h] = (t, T, I)
    rows = []
    for fp in sorted(glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')),
                     key=lambda p: parse_name(p)[1]):
        with h5py.File(fp, 'r') as hf:
            tw = (read_start_time(hf) - t0).total_seconds()
            tau, g2, g2e, q = read_g2(hf)
        noisy = np.median(g2e[(tau > 1e-4) & (tau < 1e-2), 0]) > ERR_MAX
        if tw < 0 or noisy:
            continue
        r = fit_block(tau, g2, g2e, q)
        if r is None or r['red_chi2'] >= CHI2_MAX:
            continue
        rows.append({'tw': tw, 'fit': r, 'tau': tau, 'g2': g2, 'g2e': g2e, 'q': q, 'fp': fp,
                     'f': np.array([r['per_q'][i]['f'] for i in QI]),
                     'fe': np.array([r['per_q'][i]['f_err'] for i in QI])})
    series[h] = rows




FIG_NAMES = {'H0402': 'FigureS8_Temperature_2800.pdf', 'H0407': 'FigureS9_Temperature_2825.pdf',
             'H0404': 'FigureS10_Temperature_2850.pdf', 'H0413': 'FigureS11_Temperature_2875.pdf'}
# REFS is defined above, next to the background-coefficient fit.
MAX_SHOWN = 7                         # groups drawn in (a), (b), (d)-(g)
CONVERTING = 0.85                     # f at the lowest Q below which arrest is under way
FIG_T_SIZE = (DOUBLE_COL, 5.7)
FIG_M_SIZE = (DOUBLE_COL, 2.90)   # the main-text version: top row only
MAIN_RUN = 'H0404'               # 28.50 C, the loading shown as Figure 3
MAIN_FIG = 'Figure3_Isothermal_SAXPCS.pdf'
# Two key rows below the panels, t_w above Q.  The Q key used to be a box
# inside panel (c), where it overlapped the fast fractions at early t_w; both
# keys now live under the axes and nothing is drawn over.  LEGEND_H is the
# figure fraction reserved for the pair.
LEGEND_H = 0.072
LEGEND_H_MAIN = 0.130          # the same two rows in the shorter main-text figure
LEG_LEFT = 0.075               # the left edge the axes and both key rows share
G2_YLIM = (0.98, 1.20)


def window_band(h, tw, width=590.0):
    """Median absolute band intensity over the acquisitions of one group."""
    t, T, I = traces[h]
    m = (t >= tw) & (t < tw + width)
    return float(np.median(I[m]))


def crossing(h, temp):
    """t_w at which the RTD first reached temp on the ramp (negative: before t_w = 0)."""
    t, T, I = traces[h]
    return float(t[np.argmax(T >= temp - 0.005)])


# --- the numbers the text rests on ---
print(f'\n{"run":19s} {"T":>5s} {"groups":>6s} {"chi2":>9s} {"f0 first->last":>14s} '
      f'{"ordered":>7s} {"t(f0=0.5), from t_w=0 | from 28.0 C, mid":>40s}')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    f0 = np.array([r['f'][0] for r in rows])
    tw = np.array([r['tw'] for r in rows])
    arr = [r for r in rows if r['f'][0] < CONVERTING]
    ordered = sum(all(r['f'][i] < r['f'][i + 1] for i in range(4)) for r in arr)
    below = np.nonzero(f0 < 0.5)[0]
    if len(below) and below[0] > 0:
        j = below[0]
        th = np.interp(0.5, [f0[j], f0[j - 1]], [tw[j], tw[j - 1]])
        thm = np.interp(0.5, [f0[j], f0[j - 1]], [tw[j] + 295, tw[j - 1] + 295]) - crossing(h, 28.0)
        th_txt = f'{th:8.0f} | {thm:8.0f}'
    else:
        th_txt = '   never |    never'
    chi = [r['fit']['red_chi2'] for r in rows]
    print(f'{run} {target:5.2f} {len(rows):6d} {min(chi):4.2f}-{max(chi):4.2f} '
          f'{f0[0]:5.2f}->{f0[-1]:5.2f}   {ordered:3d}/{len(arr):<3d} {th_txt:>40s}')

LATE = 0.70                           # f at the lowest Q that separates the two stages
print(f'\nD0 (1e7 A^2/s) in the early (f0 >= {LATE}) and late (f0 < {LATE}) stages:')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    dm = np.array([r['fit']['D0'] for r in rows]) / 1e7
    f0 = np.array([r['f'][0] for r in rows])
    early, late = dm[f0 >= LATE], dm[f0 < LATE]
    txt = f'  {run}: early {early[0]:.1f} -> {early[-1]:.1f} ({early[0] / early[-1]:.1f}x)'
    if len(late) > 1:
        txt += (f'; late {late.min():.1f}-{late.max():.1f} ({late.max() / late.min():.1f}x) '
                f'while 1-f0 goes {1 - f0[f0 < LATE].max():.2f} -> {1 - f0[-1]:.2f}')
    print(txt)

print('\nthe same structure, mobile or arresting: I(Q) at the first converting group over')
print('I(Q) of the last 28.00 C group, which is still mobile')
q_ref, I_ref = abs_iq(series['H0402'][-1]['fp'])
for h in ['H0407', 'H0404', 'H0413']:
    r = next(r for r in series[h] if r['f'][0] < CONVERTING)
    q, I = abs_iq(r['fp'])
    rat = [np.interp(x, q, I) / np.interp(x, q_ref, I_ref) for x in (0.004, 0.006, 0.008, 0.012, 0.02, 0.028)]
    print(f'  {RUNS[h][0]} t_w {r["tw"]:.0f} s, f0 {r["f"][0]:.2f}: ' + ' '.join(f'{x:.2f}' for x in rat))

print('\nstructure against arrest (absolute band intensity over the group, mm^-1):')
for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    conv = [r for r in rows if r['f'][0] < CONVERTING]
    if not conv:
        I_end = window_band(h, rows[-1]['tw'])
        print(f'  {run}: no arrest; band {window_band(h, 0.0):.2f} at t_w=0 -> {I_end:.2f} at {rows[-1]["tw"]:.0f} s '
              f'({I_end / window_band(h, 0.0):.0f}x), f0 >= {min(r["f"][0] for r in rows):.2f} throughout')
        continue
    a, b = conv[0], conv[-1]
    Ia, Ib = window_band(h, a['tw']), window_band(h, b['tw'])
    qa, ia = abs_iq(a['fp'])
    _, ib = abs_iq(b['fp'])
    ratio = '  '.join(f'{x:.3f}:{np.interp(x, qa, ib) / np.interp(x, qa, ia):.2f}' for x in (0.004, 0.008, 0.02))
    print(f'  {run}: {a["tw"]:.0f} -> {b["tw"]:.0f} s: 1-f0 {1 - a["f"][0]:.2f} -> {1 - b["f"][0]:.2f}, '
          f'band {Ia:.2f} -> {Ib:.2f} ({Ib / Ia:.2f}x); I(Q) ratio last/first at Q {ratio}')

print('\nramp: band intensity relative to the 20-22 C level of each loading')
for h, (run, sample, target) in RUNS.items():
    t, T, I = traces[h]
    base = np.median(I[T < 22.0])
    vals = [f'{x:.2f}C:{np.median(I[(T >= x - 0.13) & (T < x + 0.13) & (t < 0)]) / base:.1f}'
            for x in (26.0, 27.0, 27.5, 27.75, 28.0) if x < target]
    print(f'  {run}: 20-22 C {base:.4f} mm^-1; ' + '  '.join(vals))

print('\n6 C references, absolute and background-subtracted (mm^-1):')
for h in RUNS:
    q, I = abs_iq(os.path.join(DATA, f'Average_{REFS[h]}_att00_Rq0_00001_00020_results.hdf'))
    flat = (q > 0.006) & (q < 0.03)
    print(f'  {REFS[h]}: mean over 0.006-0.03 {np.mean(I[flat]):.4f} +/- {np.std(I[flat]):.4f}')


# --- FIGURES S11-S14 ---
apply_style()


def every_other_decade(x, _pos):
    power = int(round(np.log10(x)))
    return rf'$10^{{{power}}}$' if power % 2 == 0 else ''


for h, (run, sample, target) in RUNS.items():
    rows = series[h]
    shown = rows if len(rows) <= MAX_SHOWN else rows[::2]
    if shown[-1] is not rows[-1]:
        shown = shown + [rows[-1]]
    t_lo, t_hi = rows[0]['tw'], rows[-1]['tw']

    def tcol(tw):
        x = 0.5 if t_hi == t_lo else (tw - t_lo) / (t_hi - t_lo)
        return plt.cm.plasma(0.10 + 0.78 * x)

    # the first group on target, drawn black in (a) when it is too noisy to fit
    t0 = load_trace(run)[3]
    offsets = {}
    for fp in glob.glob(os.path.join(DATA, f'Average_{run}_*_results.hdf')):
        with h5py.File(fp, 'r') as hf:
            offsets[fp] = abs((read_start_time(hf) - t0).total_seconds())
    zero = min(offsets, key=offsets.get)
    if offsets[zero] > 60.0:
        zero = None
    zero_fitted = any(r['fp'] == zero for r in rows)

    # Loading MAIN_RUN is drawn twice: the seven-panel SI figure, and a
    # three-panel version that is Figure 3 of the main text.  The main-text
    # one carries only (a)-(c); the correlation functions at the four other Q
    # bins stay in the SI, where there is room for them.
    for MAIN in ([False, True] if h == MAIN_RUN else [False]):
        fig = plt.figure(figsize=FIG_M_SIZE if MAIN else FIG_T_SIZE)
        top = GridSpec(1, 3, figure=fig, left=0.075, right=0.925,
                       top=0.90 if MAIN else 0.955,
                       bottom=0.325 if MAIN else 0.595, wspace=0.42)
        ax1, ax2, ax3 = (fig.add_subplot(top[0, k]) for k in range(3))
        axq = []
        if not MAIN:
            bot = GridSpec(1, 4, figure=fig, left=0.075, right=0.925, top=0.470,
                           bottom=0.088 + LEGEND_H, wspace=0.08)
            axq = [fig.add_subplot(bot[0, k]) for k in range(4)]
            for a in axq[1:]:
                a.sharey(axq[0])
        label_panels([ax1, ax2, ax3])
        if axq:
            label_panels(axq, dx=-14, start=3)   # no y tick labels to clear there

        # (a) SAXS.  Drawn from the lowest fitted Q bin outward.  Below it the
        # points sit on the edge of the beamstop shadow, where the empty cell and
        # the sample no longer share a shape and the subtraction leaves a spike
        # that is an artefact of the subtraction rather than a measured intensity.
        key = []
        q, I = abs_iq(os.path.join(DATA, f'Average_{REFS[h]}_att00_Rq0_00001_00020_results.hdf'))
        m = (I > 0) & (q >= FIT_BAND[0])
        ax1.plot(q[m], I[m], 's', color=COLOR_6C, ms=MS, mfc='none', mew=MEW)
        lo, hi = I[m].min(), I[m].max()
        curves = [(r['fp'], tcol(r['tw']), f'{r["tw"]:.0f}') for r in shown]
        if zero is not None and not zero_fitted:
            curves.insert(0, (zero, 'black', '0'))
        for fp, c, label in curves:
            q, I = abs_iq(fp)
            m = (I > 0) & (q >= FIT_BAND[0])
            ax1.plot(q[m], I[m], 'o', color=c, ms=MS, mfc='none', mew=MEW)
            lo, hi = min(lo, I[m].min()), max(hi, I[m].max())
            key.append(Line2D([], [], marker='o', ls='none', color=c, ms=MS, mfc='none', mew=MEW, label=label))
        ax1.set_xscale('log')
        ax1.set_yscale('log')
        ax1.set_xlim(3.3e-3, 3.6e-2)
        q_log_ticks(ax1)
        ax1.set_ylim(lo * 0.45, hi * 3.0)
        q6, I6 = abs_iq(os.path.join(DATA, f'Average_{REFS[h]}_att00_Rq0_00001_00020_results.hdf'))
        plateau = np.median(I6[(q6 > 0.006) & (q6 < 0.03)])
        ax1.text(0.012, plateau * 2.4, '6 °C ref', ha='center', va='bottom', color=COLOR_6C)
        ax1.set_xlabel(r'$Q$ ($\AA^{-1}$)')
        ax1.set_ylabel(r'$I(Q)$ (mm$^{-1}$)')
        add_minor_grid(ax1)

        # (b) and (d)-(g): g2 with fits
        for ax, qi in [(ax2, 0)] + list(zip(axq, [1, 2, 3, 4])):
            for r in shown:
                c = tcol(r['tw'])
                m = r['tau'] > 0
                ax.errorbar(r['tau'][m], r['g2'][m, qi], yerr=r['g2e'][m, qi], fmt='o', color=c,
                            mfc='none', ms=MS - 0.5, mew=LW_THIN, capsize=1.0, elinewidth=LW_THIN,
                            capthick=LW_THIN, alpha=0.75, zorder=2)
                fr, pq = r['fit'], r['fit']['per_q'][qi]
                tt = np.logspace(np.log10(r['tau'][m].min()), np.log10(r['tau'][m].max()), 200)
                ax.plot(tt, xpcs_fit.double_exp(tt, pq['tau_fast'], pq['f'], pq['tau_slow'],
                                                fr['p1'], fr['p2']), color=c, lw=LW_DATA, zorder=3)
            # The contrast the fit uses, drawn so the reader can see how close each
            # correlation function comes to full arrest without reading a number
            # off the caption.  It is 1 + BETA, the plateau the static reference
            # reaches between 0.0085 and 0.0165 A^-1, NOT that reference's level
            # at this Q: the standard is scratched, so its apparent contrast in
            # these bins is depressed by up to 12% and a line drawn there would
            # sit below data that are perfectly ordinary (Figure S6).
            ax.axhline(1.0 + BETA, color='0.35', ls='--', lw=LW_DATA * 1.7, zorder=1)
            ax.set_xscale('log')
            ax.set_xlim(9e-6, 5.0)
            ax.set_ylim(*G2_YLIM)
            ax.set_yticks([1.00, 1.05, 1.10, 1.15])
            ax.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
            ax.xaxis.set_major_formatter(FuncFormatter(every_other_decade))
            ax.set_xlabel(r'Delay Time, $\Delta t$ (s)')
            add_minor_grid(ax)
            ax.text(0.96, 0.95, rf'$Q = {rows[0]["q"][qi]:.5f}\ \AA^{{-1}}$', transform=ax.transAxes,
                    ha='right', va='top')
        ax2.annotate(r'$1+\beta$', xy=(0.985, 1.0 + BETA),
                     xycoords=ax2.get_yaxis_transform(), xytext=(0, 4),
                     textcoords='offset points', ha='right', va='bottom', color='0.35')
        ax2.set_ylabel(r'$g_2$')
        if axq:
            axq[0].set_ylabel(r'$g_2$')
            for a in axq[1:]:
                plt.setp(a.get_yticklabels(), visible=False)

        # (c) f against t_w, colour = t_w, marker = Q bin; D0 on the right axis
        for i, mk in zip(QI, Q_MARKERS):
            x = [r['tw'] for r in rows]
            ax3.plot(x, [r['f'][i] for r in rows], '-', color='0.75', lw=LW_THIN, zorder=1)
            for r in rows:
                ax3.errorbar(r['tw'], r['f'][i], yerr=r['fe'][i], marker=mk, color=tcol(r['tw']),
                             ms=MS_SPARSE - 1, capsize=1.0, mfc='none', mew=MEW,
                             elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
        f_lo = min(min(r['f'] - r['fe']) for r in rows)
        f_hi = max(max(r['f'] + r['fe']) for r in rows)
        span = max(f_hi - f_lo, 0.2)
        ax3.set_ylim(f_lo - 0.08 * span, f_hi + 0.42 * span)
        step = 0.2 if span > 0.5 else 0.1
        ax3.set_yticks([v for v in np.arange(0.0, 1.0001, step) if f_lo - 0.08 * span <= v])
        pad = 0.12 * max(t_hi - t_lo, 600)
        ax3.set_xlim(t_lo - pad, t_hi + pad)
        ax3.set_xlabel(r'$t_w$ (s)')
        ax3.set_ylabel('Fast Fraction, $f$')
        add_minor_grid(ax3)
        qh = [Line2D([], [], marker=mk, ls='none', mfc='none', mec='k', mew=MEW,
                     ms=MS_SPARSE - 1, label=f'{1e3 * rows[0]["q"][i]:.2f}')
              for i, mk in zip(QI, Q_MARKERS)]
        ax3r = ax3.twinx()
        yd = [r['fit']['D0'] / 1e7 for r in rows]
        ed = [r['fit']['D0_err'] / 1e7 for r in rows]
        ax3r.errorbar([r['tw'] for r in rows], yd, yerr=ed, marker='*', ls='none', color='k',
                      ms=MS_SPARSE + 2.0, capsize=1.0, mfc='none', mew=MEW, elinewidth=LW_THIN,
                      capthick=LW_THIN, zorder=4)
        ax3r.set_ylim(0, max(y + e for y, e in zip(yd, ed)) * 1.30)
        ax3r.set_ylabel(r'$D_0$ ($10^{7}\ \AA^{2}\,$s$^{-1}$)')

        # Two rows, each one line: the label sits inline ahead of its entries,
        # carried by a handle that draws nothing.  Both rows are anchored at the
        # same left edge as the axes, so they align with each other.
        row_h = LEGEND_H_MAIN if MAIN else LEGEND_H
        for handles, title, yy in ((key, r'$t_w$ (s)', row_h * 0.52),
                                   (qh, r'$Q$ ($10^{-3}\ \AA^{-1}$)', 0.004)):
            hh = [Line2D([], [], ls='none', marker='', label=title)] + list(handles)
            lg = fig.legend(handles=hh, loc='lower left',
                            bbox_to_anchor=(LEG_LEFT, yy), ncol=len(hh),
                            handlelength=1.0, handletextpad=0.35, columnspacing=1.0,
                            borderaxespad=0.0, borderpad=0.3, labelspacing=0.2)
            lg.get_frame().set_linewidth(LW_THIN)
        save_fig(fig, MAIN_FIG if MAIN else FIG_NAMES[h])

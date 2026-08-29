"""Combined SAXS + XPCS analysis, reading the averaged HDF files in data/.

This one script writes three of the paper's figures: Figure 3 of the main text,
Figure S9 and Figure S3.

Figure 3, three panels:
  1. SAXS I(q) for B0146, the first B0147 file (frames 1-200), and the last five
     B0147 files (frames 801-1313), with the D0138 buffer subtracted.
  2. XPCS g2(tau) at one q for the first + last-five B0147 files, with fits.
  3. Fitted fast fraction f vs elapsed time, per q bin.
Figure S9 (2x2 panels) plots the fit parameters: (a) shared exponents p1, p2 vs
elapsed time, (b) tau_fast vs Q, (c) tau_slow vs Q, and (d) the power-law
scaling exponents gamma_fast, gamma_slow obtained by fitting each elapsed
time's tau(Q) in (b)/(c) to tau = A * Q**gamma.  Figure S3 (2 panels) documents
the absolute-cross-section calibration: the ion-chamber -> photon linear fit
and the air transmission (see the ABSOLUTE SCATTERING CROSS-SECTION section
below and abs_xsec.py).

The SAXS panel is still put on an absolute scale (d(Sigma)/d(Omega)) via a
coefficient computed PER FILE from that file's own range-averaged ion-chamber
readings, so a drift in incident flux across the time series is handled
correctly.  The axis carries the units (cm^-1); the calibration behind them is
in abs_xsec_coef() and in the long note in the ABSOLUTE SCATTERING CROSS-SECTION
section.

Colour encodes elapsed time and is CONSISTENT across all three panels (e.g. the
7863 s dataset is the same gold everywhere).  In the fit panel the marker SHAPE
encodes the q bin.  B0146 is a 6 C reference taken before the 30 C isothermal
run; its acquisition time is not a time origin, so it appears in the SAXS panel
only.  Because colour means the same thing everywhere, ONE elapsed-time key
serves all three panels, and it therefore belongs to the figure rather than to
any one panel: it is a single-row box along the bottom.  Keeping it out of the
axes means no panel has to carry headroom for it, so all three y ranges sit
close to their data.  The 6 C reference is NOT one of the elapsed times and
appears in panel (a) only, so it is keyed separately in that panel's
bottom-left corner, beside its own curve.  The marker-shape -> q key stays
inside panel (c), flattened to three columns so it too costs little headroom.

Figure geometry follows the ACS figure-preparation guidelines through
common/acs_style.py: all three multi-panel figures are double-column (7.0 in = the
504 pt ACS maximum) and every character is 8 pt Arial.  Axis limits are set a
clear margin outside the data everywhere, so no point is drawn on top of a
frame.

Fit model and shared exponents
------------------------------
The g2 model is a double stretched-exponential (Siegert form):

    g2 = contrast * ( f e^-(tau/tau_fast)^p1 + (1-f) e^-(tau/tau_slow)^p2 )^2 + 1

with the contrast fixed at beta = 0.13042, the instrumental value measured on a
static reference by contrast_calibration.py (Figure S4), and the baseline fixed
at 1.  For each
elapsed time all fitted q bins are fit SIMULTANEOUSLY (a global fit): the
stretching exponents p1 (fast) and p2 (slow) are SHARED across q -- they depend
only on elapsed time -- while tau_fast, f and tau_slow are independent per q.
This is the physically motivated constraint that the relaxation-shape exponents
are a property of the sample state at a given age, not of the q bin, and it
stabilises the otherwise poorly-conditioned slow mode (tau_slow lies beyond the
~1.5 s delay window, so its shape cannot be pinned q-by-q).

Uncertainties
-------------
The fit minimises error-weighted residuals (g2_model - g2)/g2_err using the
g2_err stored in the averaged files directly (absolute_sigma convention).  The
parameter covariance is inv(J^T J) at the solution; 1-sigma errors on f, p1 and
p2 are the square roots of its diagonal.  Because p1/p2 are shared, their
uncertainty is correctly propagated into f (f errors are larger than a per-q
fixed-exponent fit would report -- that difference is the exponent systematic,
now folded in honestly).
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
from matplotlib.lines import Line2D
from matplotlib.ticker import (FixedLocator, FixedFormatter, NullFormatter,
                               FuncFormatter, LogLocator)
from scipy.optimize import least_squares

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_SPARSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)

# --- PARAMETERS ---
data_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "data")

# --- ACS FIGURE GEOMETRY (see common/acs_style.py) ---
# All three have two or more panels, so all three are double-column figures at
# the 7.0 in (504 pt) ACS maximum.  Depths are far inside the 9.167 in cap.
FIG1_SIZE = (DOUBLE_COL, 3.15)        # 3 panels + the shared key row beneath them
FIG1_LEGEND_H = 0.075                # fraction of the height reserved at the
                                     # bottom for that shared key
FIG2_SIZE = (DOUBLE_COL, 4.2)        # 2x2 fit-parameter panels
FIG3_SIZE = (DOUBLE_COL, 2.9)        # 2 calibration panels
# Panel (a) carries no key of its own except the one-line 6 C reference in the
# bottom-left corner, so its limits sit just clear of the data: a little over a
# third of a decade below, half a decade above.
SAXS_YLO_PAD, SAXS_YHI_PAD = 0.45, 3.0

# --- DATASET SELECTION ---
BACKGROUND_HEADER = 'D0138'          # buffer, subtracted as background in the SAXS panel
XPCS_HEADER = 'B0147'                # header used for the g2 / fit panels
# The absolute cross-section coefficients coef_sam / coef_buf are no longer
# hard-coded: they are computed per file from that file's own (range-averaged)
# ion-chamber readings by abs_xsec_coef(), defined in the ABSOLUTE SCATTERING
# CROSS-SECTION section below.
N_LAST = 5                           # number of last (highest-frame) files
target_q_idx = 0                     # q bin shown in the g2 panel
fit_q_indices = [0, 1, 2, 3, 4]      # q bins fitted for the f-vs-time panel
Q_MARKERS = ['o', 's', '^', 'D', 'v']   # one marker per fit q bin
CHI2_MAX = 10.0                      # skip fits worse than this reduced chi^2.
                                     # A guard only: the five fitted files give
                                     # 0.78-1.89, so it never fires.  The liquid
                                     # t=0 file, which it was written for, is not
                                     # in the fit set at all (chi^2 ~ 2e3 if it
                                     # were: it does not obey the arrested model)
PHI_AVERAGE = True                   # azimuthally average I(q) over phi sectors
# XPCS elapsed-time colours: the N_LAST times are mapped to evenly-spaced
# positions of this colormap (maximises separation so times are easy to tell
# apart), independent of the SAXS-only SI figure.  The 0 s and 6 C reference
# datasets get their own distinct colours (below), outside this scale.
XPCS_CMAP = plt.cm.plasma
CMAP_LO, CMAP_HI = 0.10, 0.88        # colormap span used for the N_LAST times
COLOR_0S = 'black'                   # first B0147 file (elapsed = 0 s)
COLOR_6C = '#1f77b4'                 # B0146 (6 C reference, before isothermal)

# --- FIT MODEL ---
# The model, the measured contrast and the global fit live in xpcs_fit.py so
# that Figure 3b, Figure S9 and Figure S10 provably share one implementation.
from xpcs_fit import (CONTRAST, BASELINE, double_exp, fit_g2_global,  # noqa: E402
                      PQ_P0, PQ_LO, PQ_HI, P_EXP_P0, P_EXP_LO, P_EXP_HI)
contrast = CONTRAST                  # local alias used in the panels below


# ============================================================================
# ABSOLUTE SCATTERING CROSS-SECTION
# ============================================================================
# The calibration constants, the IC->photon linear fit, the air transmission
# and the per-file coefficient function abs_xsec_coef() all live in the shared
# module abs_xsec.py.  saxpcs.py and saxs_evolution.py both import from it so
# the absolute scale is computed identically in every figure.  See abs_xsec.py
# for the full derivation, its provenance and the unit notes.
from abs_xsec import (                                     # noqa: E402
    CAL_A, CAL_B, CAL_UPIC, CAL_PHOTONS, CAL_CROP,
    AIR_TRANS_SERIES, AIR_TRANSMISSION, AIR_TRANS_STD,
    INCIDENT_PATH, TRANSMITTED_PATH, INV_MM_TO_INV_CM,
    abs_xsec_coef, calibration_summary,
)

print(calibration_summary())


# --- HDF field locations ---
START_TIME_PATH = '/entry/start_time'
TIME_FORMAT     = '%Y-%m-%d %H:%M:%S'
FRAME_TIME_PATH = '/entry/instrument/detector_1/frame_time'
DELAY_PATH      = '/xpcs/multitau/delay_list'
G2_PATH         = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH     = '/xpcs/multitau/normalized_g2_err'
DYN_Q_PATH      = '/xpcs/qmap/dynamic_v_list_dim0'
SAXS_PATH       = '/xpcs/temporal_mean/scattering_1d'
STATIC_MAP_PATH = '/xpcs/qmap/static_index_mapping'
STATIC_Q_PATH   = '/xpcs/qmap/static_v_list_dim0'
STATIC_PHI_PATH = '/xpcs/qmap/static_v_list_dim1'

_name_re = re.compile(r'Average_([A-Za-z]\d+)_.*?_(\d+)_(\d+)_results')


# --- READ HELPERS ---
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


def read_g2(hf):
    t0 = hf[FRAME_TIME_PATH][()]
    t0 = t0.item() if isinstance(t0, np.ndarray) else t0
    tau = hf[DELAY_PATH][()] * t0
    tau = tau[:, 0] if tau.ndim > 1 else tau
    return tau, hf[G2_PATH][()], hf[G2_ERR_PATH][()], hf[DYN_Q_PATH][()]


def fit_powerlaw(Q, tau, tau_err):
    """Weighted log-log fit of tau = A * Q**gamma at one elapsed time.

    Minimises the error-weighted residual of ln(tau) vs ln(Q) (a weighted
    linear regression), with sigma_ln(tau) = tau_err / tau (propagated from
    the tau uncertainty). Returns (gamma, gamma_err, A), with gamma_err from
    the weighted-least-squares covariance (absolute_sigma convention) and A
    the fitted prefactor (for overlaying the fit line).
    """
    x = np.log(np.asarray(Q, dtype=float))
    y = np.log(np.asarray(tau, dtype=float))
    sigma_y = np.asarray(tau_err, dtype=float) / np.asarray(tau, dtype=float)
    w = 1.0 / sigma_y**2
    X = np.vstack([x, np.ones_like(x)]).T
    cov = np.linalg.inv(X.T @ (w[:, None] * X))
    gamma, intercept = cov @ (X.T @ (w * y))
    gamma_err = np.sqrt(cov[0, 0])
    return gamma, gamma_err, np.exp(intercept)


# --- DISCOVER FILES ---
# every range average in data/, including the thermal-cycle groups that belong to
# Figure S6; the headers wanted here are selected out of by_header below
file_paths = sorted(glob.glob(os.path.join(data_dir, 'Average_*.hdf')))
assert file_paths, f'no HDF files found in {data_dir}'
by_header, start_times = {}, {}
for fp in file_paths:
    header = parse_name(fp)[0]
    with h5py.File(fp, 'r') as hf:
        start_times[fp] = read_start_time(hf)
    by_header.setdefault(header, []).append(fp)
for h in by_header:
    by_header[h].sort(key=lambda p: parse_name(p)[1])

b0147 = by_header[XPCS_HEADER]
first_file = b0147[0]                          # first B0147 file (frames 1-200)
xpcs_files = b0147[-N_LAST:]                    # g2 / fit: last N only (NOT 1-200)
saxs_files = [first_file] + xpcs_files          # SAXS panel: first + last N
t_ref = start_times[first_file]                 # time origin = first B0147 file


def elapsed(fp):
    return (start_times[fp] - t_ref).total_seconds()


# XPCS elapsed-time colours: map the N_LAST files (sorted by time) to evenly
# spaced colormap positions so consecutive times are maximally distinguishable.
# The same mapping is reused in every XPCS panel (and in g2_grid_SI.py).
_xpcs_sorted = sorted(xpcs_files, key=elapsed)
_xpcs_pos = np.linspace(CMAP_LO, CMAP_HI, len(_xpcs_sorted))
_xpcs_color = {fp: XPCS_CMAP(p) for fp, p in zip(_xpcs_sorted, _xpcs_pos)}


def ecolor(fp):
    """Colour for an XPCS elapsed time; 0 s and the 6 C ref use distinct colours."""
    if fp == first_file:
        return COLOR_0S
    return _xpcs_color[fp]


# elapsed-time value -> colour (xpcs_files only; used where fp isn't in scope)
time_color = {elapsed(fp): ecolor(fp) for fp in xpcs_files}


print('SAXS files:', [f'{parse_name(f)[1]}-{parse_name(f)[2]} ({elapsed(f):.0f} s)' for f in saxs_files])
print('XPCS files:', [f'{parse_name(f)[1]}-{parse_name(f)[2]} ({elapsed(f):.0f} s)' for f in xpcs_files])

# --- GLOBAL FITS (one per elapsed time; p1, p2 shared across q) ---
# Fit each XPCS file once and reuse the result in every panel.
g2_data = {}   # fp -> (tau, g2, g2_err, q_vals)
fits = {}      # fp -> result dict from fit_g2_global
for fp in xpcs_files:
    with h5py.File(fp, 'r') as hf:
        g2_data[fp] = read_g2(hf)
    tau, g2, g2_err, q_vals = g2_data[fp]
    fits[fp] = fit_g2_global(tau, g2, g2_err, fit_q_indices)
    r = fits[fp]
    if r is not None:
        fmt = ' '.join(f'f{qi}={r["per_q"][qi]["f"]:.3f}' for qi in r['per_q'])
        print(f'  fit {elapsed(fp):5.0f} s: p1={r["p1"]:.3f}+/-{r["p1_err"]:.3f} '
              f'p2={r["p2"]:.3f}+/-{r["p2_err"]:.3f} chi2={r["red_chi2"]:.2f}  {fmt}')

# --- FIGURE ---
apply_style()
fig, (ax1, ax2, ax3) = plt.subplots(1, 3, figsize=FIG1_SIZE)
label_panels((ax1, ax2, ax3))

# ============================================================
# PANEL 1: SAXS I(q), background-subtracted (B0146 + first + last N B0147)
# ============================================================
# Background buffer curve AND its own absolute coefficient (coef_buf).  Both are
# taken from the D0138 buffer file; coef_buf multiplies the buffer term in every
# subtraction below.
bg_I = None
coef_buf = None
if BACKGROUND_HEADER in by_header:
    with h5py.File(by_header[BACKGROUND_HEADER][0], 'r') as hf:
        _, bg_I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_buf = abs_xsec_coef(hf)
    print(f'coef_buf ({BACKGROUND_HEADER}) = {coef_buf:.3e}')

time_handles = []                         # the six elapsed times -> the key row
ref_handles = []                          # the 6 C reference -> its own small key
saxs_I_lo, saxs_I_hi = np.inf, 0.0        # data range, for the panel-(a) headroom
for fp in by_header.get('B0146', []):                        # 6 C reference (distinct)
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                          # this file's own coefficient
    print(f'coef_sam (B0146 {parse_name(fp)[1]}-{parse_name(fp)[2]}) = {coef_sam:.3e}')
    I = coef_sam * I - coef_buf * bg_I if bg_I is not None else coef_sam * I
    I = INV_MM_TO_INV_CM * I          # abs_xsec_coef() is mm^-1; the axis is cm^-1
    pos = I > 0
    saxs_I_lo, saxs_I_hi = min(saxs_I_lo, I[pos].min()), max(saxs_I_hi, I[pos].max())
    ax1.plot(q[pos], I[pos], color=COLOR_6C, marker='s', ls='none', ms=MS,
             mfc='none', mew=MEW)
    ref_handles.append(Line2D([], [], color=COLOR_6C, marker='s', ls='none',
                              mfc='none', mew=MEW, ms=MS,
                              label='6 °C ref'))

for fp in saxs_files:                                        # first (0 s) + last N B0147
    with h5py.File(fp, 'r') as hf:
        q, I = read_saxs_iq(hf, PHI_AVERAGE)
        coef_sam = abs_xsec_coef(hf)                          # this file's own coefficient
    print(f'coef_sam (B0147 {parse_name(fp)[1]}-{parse_name(fp)[2]}) = {coef_sam:.3e}')
    I = coef_sam * I - coef_buf * bg_I if bg_I is not None else coef_sam * I
    I = INV_MM_TO_INV_CM * I          # abs_xsec_coef() is mm^-1; the axis is cm^-1
    pos = I > 0
    saxs_I_lo, saxs_I_hi = min(saxs_I_lo, I[pos].min()), max(saxs_I_hi, I[pos].max())
    color = ecolor(fp)
    ax1.plot(q[pos], I[pos], color=color, marker='o', ls='none', ms=MS,
             mfc='none', mew=MEW)
    time_handles.append(Line2D([], [], color=color, marker='o', ls='none',
                               mfc='none', mew=MEW, ms=MS,
                               label=f'{elapsed(fp):.0f} s'))

ax1.set_xscale('log')
ax1.set_yscale('log')
ax1.set_xlabel(r'$Q$ ($\AA^{-1}$)')
ax1.set_ylabel(r'$I(Q)$ (cm$^{-1}$)')
add_minor_grid(ax1)
# Q spans only ~1 decade (0.0033-0.034), so the automatic log locator labels a
# single tick (10^-2).  Label three well-separated decimal positions instead --
# far enough apart that the labels cannot touch at 8 pt.  The limits sit well
# outside the data so no point sits on the frame.
ax1.set_xlim(2.6e-3, 4.4e-2)
ax1.xaxis.set_major_locator(FixedLocator([4e-3, 1e-2, 3e-2]))
ax1.xaxis.set_major_formatter(FixedFormatter(['0.004', '0.01', '0.03']))
# Minor ticks only at integer multiples of the decade (2, 3, ... x 10^n).
# Half-decade positions such as 0.015 or 0.025 put grid lines at values a reader
# cannot name, which makes the grid harder to read rather than easier.
ax1.xaxis.set_minor_locator(FixedLocator([3e-3, 5e-3, 6e-3, 7e-3, 8e-3, 9e-3,
                                          2e-2, 4e-2]))
ax1.xaxis.set_minor_formatter(NullFormatter())

# The elapsed-time key lives at the bottom of the FIGURE (below), so panel (a)
# needs no headroom for it.  Only the one-line 6 C reference key sits inside,
# in the bottom-left corner next to the curve it labels: that corner is empty
# because I(Q) falls to the right.
ax1.set_ylim(saxs_I_lo * SAXS_YLO_PAD, saxs_I_hi * SAXS_YHI_PAD)
leg_r = ax1.legend(handles=ref_handles, loc='lower left', borderaxespad=0.3,
                   handlelength=1.2, handletextpad=0.3)
leg_r.get_frame().set_linewidth(LW_THIN)

# ============================================================
# PANEL 2: g2(tau) at target q, with fits
# ============================================================
q_val_target = None
for fp in xpcs_files:
    tau, g2, g2_err, q_vals = g2_data[fp]
    q_val_target = q_vals[target_q_idx]
    color = ecolor(fp)
    m = tau > 0
    # Data drawn at the THIN end of the allowed range (0.5 pt = 0.18 mm edges)
    # and below the fits in z, so the fit line reads as a line and not as one
    # more marker outline of the same colour.
    ax2.errorbar(tau[m], g2[m, target_q_idx], yerr=g2_err[m, target_q_idx], fmt='o',
                 color=color, mfc='none', markersize=MS, mew=LW_THIN,
                 capsize=1.5, elinewidth=LW_THIN, capthick=LW_THIN,
                 alpha=0.75, linestyle='none', zorder=2)
    r = fits[fp]
    if r is not None and target_q_idx in r['per_q'] and r['red_chi2'] < CHI2_MAX:
        pq = r['per_q'][target_q_idx]
        tt = tau[m]
        t_fit = np.logspace(np.log10(tt.min()), np.log10(tt.max()), 200)
        g2_fit = double_exp(t_fit, pq['tau_fast'], pq['f'], pq['tau_slow'],
                            r['p1'], r['p2'])
        # LW_DATA (1.0 pt) against 0.6 pt marker outlines is enough contrast on
        # its own -- deliberately NO white casing under the line, which would
        # blank out the very data points the fit is meant to be judged against.
        ax2.plot(t_fit, g2_fit, color=color, lw=LW_DATA, solid_capstyle='round',
                 zorder=3)

ax2.set_xscale('log')
ax2.set_xlabel(r'Delay Time, $\tau$ (s)')
ax2.set_ylabel(r'$g_2$')
ax2.set_ylim(1.0, 1.18)
# two-decimal ticks: '1.000'-style labels are wide enough at 8 pt to squeeze
# the three panels together, and the extra digit carries no information here.
ax2.yaxis.set_major_locator(FixedLocator([1.00, 1.05, 1.10, 1.15]))
ax2.yaxis.set_major_formatter(FixedFormatter(['1.00', '1.05', '1.10', '1.15']))
add_minor_grid(ax2)
# tau covers ~5 decades; label every other one so the exponents stay legible
# and well separated, keeping a tick (and grid line) on every decade.  The
# limits leave a clear margin either side of the shortest / longest delay.
ax2.set_xlim(9e-6, 5.0)
ax2.xaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
ax2.xaxis.set_major_formatter(FuncFormatter(
    lambda x, _pos: rf'$10^{{{int(round(np.log10(x)))}}}$'
                    if int(round(np.log10(x))) % 2 == 0 else ''))

# ============================================================
# PANEL 3: fitted fast fraction f vs elapsed time (colour=time, marker=q)
# ============================================================
fit_rows = []
q_val_of = {}                                    # q_idx -> Q value (for legends)
for fp in xpcs_files:
    r = fits[fp]
    if r is None or r['red_chi2'] >= CHI2_MAX:   # drop unreliable global fits
        continue
    _, _, _, q_vals = g2_data[fp]
    for q_idx, pq in r['per_q'].items():
        q_val_of[q_idx] = q_vals[q_idx]
        fit_rows.append({'q_index': q_idx, 'q_val': q_vals[q_idx],
                         'elapsed': elapsed(fp), 'f': pq['f'], 'f_err': pq['f_err']})

for q_idx, marker in zip(fit_q_indices, Q_MARKERS):
    rows = sorted([r for r in fit_rows if r['q_index'] == q_idx], key=lambda r: r['elapsed'])
    if not rows:
        continue
    xs = [r['elapsed'] for r in rows]
    ys = [r['f'] for r in rows]
    es = [r['f_err'] for r in rows]
    ax3.plot(xs, ys, '-', color='0.75', lw=LW_THIN, zorder=1)        # per-q guide line
    for x, y, e in zip(xs, ys, es):
        ax3.errorbar(x, y, yerr=e, marker=marker, color=time_color[x],
                     markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                     elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)

# Headroom above the highest fitted point, for the q key.  Three columns make
# that key two rows deep instead of five, so 0.18 of the span is enough where a
# single column needed 0.42 -- the data now fills most of the panel.
_f_hi = max(r['f'] + r['f_err'] for r in fit_rows)
_f_lo = min(r['f'] - r['f_err'] for r in fit_rows)
_span = _f_hi - _f_lo
ax3.set_ylim(_f_lo - 0.08 * _span, _f_hi + 0.18 * _span)
# widen the elapsed-time axis so the first and last groups are not on the frame
_t_lo = min(r['elapsed'] for r in fit_rows)
_t_hi = max(r['elapsed'] for r in fit_rows)
_t_pad = 0.16 * (_t_hi - _t_lo)
ax3.set_xlim(_t_lo - _t_pad, _t_hi + _t_pad)

# legend: marker shape -> q bin (open black markers).  The Q values are carried
# by a legend TITLE with the 10^-3 factor pulled out, so each entry is a short
# mantissa instead of '$Q = 0.00376 \AA^{-1}$'.  Three columns keep the box
# shallow, so it clears the descending 5039 s points on little headroom.
q_handles = [Line2D([], [], marker=mk, ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=f'{1e3 * q_val_of[qi]:.2f}')
             for qi, mk in zip(fit_q_indices, Q_MARKERS) if qi in q_val_of]
leg3 = ax3.legend(handles=q_handles, loc='upper right', ncol=3,
                  borderaxespad=0.3, title=r'$Q$ ($10^{-3}\ \AA^{-1}$)',
                  handlelength=0.8, handletextpad=0.25, columnspacing=0.5,
                  labelspacing=0.2, borderpad=0.3)
leg3.get_frame().set_linewidth(LW_THIN)

ax3.set_xlabel('Elapsed Time (s)')
ax3.set_ylabel('Fast Fraction, $f$')
add_minor_grid(ax3)

# --- shared elapsed-time key, one row along the bottom of the figure ---
# Colour carries the same meaning in all three panels, so this key belongs to
# the figure, not to any single panel.  Out here it costs no panel any y range.
leg_fig = fig.legend(handles=time_handles, loc='lower center',
                     bbox_to_anchor=(0.5, 0.0), ncol=len(time_handles),
                     handlelength=1.2, handletextpad=0.4, columnspacing=1.2,
                     borderaxespad=0.0, borderpad=0.35)
leg_fig.get_frame().set_linewidth(LW_THIN)

# Reserve the bottom strip for that key; w_pad keeps the y labels of panels (b)
# and (c) clear of the neighbouring panel's tick labels.
fig.tight_layout(pad=0.4, w_pad=1.1, rect=(0, FIG1_LEGEND_H, 1, 1))
save_fig(fig, 'Figure3_Isothermal_SAXPCS.pdf')

# ============================================================
# FIGURE S9 (2x2): (a) p1, p2 vs elapsed time;      (b) tau_fast vs Q
#                  (d) gamma_fast, gamma_slow vs t;  (c) tau_slow vs Q
# Column 0 (a, d) shares the elapsed-time x-axis; column 1 (b, c) shares the Q
# x-axis -- so only the bottom row needs x tick labels / an x-axis label.
# All XPCS colours use the same elapsed-time scale as Figure 3.
# ============================================================
fig2, ((axp, axf), (axg, axs)) = plt.subplots(2, 2, figsize=FIG2_SIZE,
                                              sharex='col')
label_panels((axp, axf, axs, axg))

# --- (a) shared stretching exponents vs elapsed time (colour = time) ---
exp_rows = [{'elapsed': elapsed(fp), 'fp': fp,
             'p1': fits[fp]['p1'], 'p1_err': fits[fp]['p1_err'],
             'p2': fits[fp]['p2'], 'p2_err': fits[fp]['p2_err']}
            for fp in xpcs_files
            if fits[fp] is not None and fits[fp]['red_chi2'] < CHI2_MAX]
exp_rows.sort(key=lambda r: r['elapsed'])
xe = [r['elapsed'] for r in exp_rows]
axp.plot(xe, [r['p1'] for r in exp_rows], '-', color='0.75', lw=LW_THIN, zorder=1)
axp.plot(xe, [r['p2'] for r in exp_rows], '-', color='0.75', lw=LW_THIN, zorder=1)
for r in exp_rows:                                   # p1 = circle, p2 = square
    axp.errorbar(r['elapsed'], r['p1'], yerr=r['p1_err'], marker='o', color=ecolor(r['fp']),
                 markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
    axp.errorbar(r['elapsed'], r['p2'], yerr=r['p2_err'], marker='s', color=ecolor(r['fp']),
                 markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
axp.axhline(1.0, color='0.6', ls=':', lw=LW_THIN)          # simple-exponential reference
p_handles = [Line2D([], [], marker='o', ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=r'$p_{\mathrm{fast}}$'),
             Line2D([], [], marker='s', ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=r'$p_{\mathrm{slow}}$')]
axp.legend(handles=p_handles, loc='upper left')
# headroom above BOTH the data and the p = 1 reference line, so the upper-left
# legend box has clear space and does not sit on the dotted line
_p_lo = min(min(r['p1'] - r['p1_err'], r['p2'] - r['p2_err']) for r in exp_rows)
_p_hi = max(max(r['p1'] + r['p1_err'], r['p2'] + r['p2_err']) for r in exp_rows)
_p_span = max(_p_hi, 1.0) - _p_lo
axp.set_ylim(_p_lo - 0.10 * _p_span, max(_p_hi, 1.0) + 0.18 * _p_span)
axp.set_ylabel('Stretching exponent')
axp.tick_params(labelbottom=False)
add_minor_grid(axp)

# --- (b, c) relaxation times vs Q, one curve per elapsed time; also fit each
# elapsed time's tau(Q) to a power law tau = A * Q**gamma for panel (d). ---
gamma_rows = []
for fp in xpcs_files:
    r = fits[fp]
    if r is None or r['red_chi2'] >= CHI2_MAX:
        continue
    _, _, _, q_vals = g2_data[fp]
    qs = sorted(r['per_q'])
    Q      = np.array([q_vals[qi] for qi in qs])
    tf     = np.array([r['per_q'][qi]['tau_fast'] for qi in qs])
    tf_err = np.array([r['per_q'][qi]['tau_fast_err'] for qi in qs])
    ts     = np.array([r['per_q'][qi]['tau_slow'] for qi in qs])
    ts_err = np.array([r['per_q'][qi]['tau_slow_err'] for qi in qs])
    f_val  = np.array([r['per_q'][qi]['f'] for qi in qs])
    f_err  = np.array([r['per_q'][qi]['f_err'] for qi in qs])
    # tau_fast only means something where a fast mode is actually detected.  At
    # the last elapsed time the lowest q bin fits f = 0.000 +/- 0.003, i.e. no
    # fast mode at all, and tau_fast there rails to a meaningless ~5 s.  Keeping
    # that point would drag the tau_fast(Q) power law to gamma = -15.  Bins whose
    # fast amplitude is under 3 sigma are therefore dropped from panel (b) and
    # from the gamma_fast fit; the slow mode carries every bin.
    det = f_val > 3 * f_err
    color  = ecolor(fp)
    lbl    = f'{elapsed(fp):.0f} s'
    axf.errorbar(Q[det], tf[det], yerr=tf_err[det], marker='o', ls='none', color=color,
                 markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, label=lbl, zorder=2)
    axs.errorbar(Q, ts, yerr=ts_err, marker='s', ls='none', color=color,
                 markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, label=lbl, zorder=2)
    if det.sum() < len(Q):
        print(f'    (t_w={elapsed(fp):.0f} s: {int((~det).sum())} q bin(s) dropped '
              f'from the fast-mode panel, f < 3 sigma)')

    if len(Q) >= 3:
        gf, gf_err, gf_A = fit_powerlaw(Q[det], tf[det], tf_err[det]) if det.sum() >= 3 \
            else (np.nan, np.nan, np.nan)
        gs, gs_err, gs_A = fit_powerlaw(Q, ts, ts_err)
        Q_line = np.linspace(Q.min(), Q.max(), 50)
        if np.isfinite(gf):
            axf.plot(np.linspace(Q[det].min(), Q[det].max(), 50),
                     gf_A * np.linspace(Q[det].min(), Q[det].max(), 50)**gf,
                     '--', color=color, lw=LW_DATA, zorder=1)
        axs.plot(Q_line, gs_A * Q_line**gs, '--', color=color, lw=LW_DATA, zorder=1)
        gamma_rows.append({'elapsed': elapsed(fp), 'fp': fp,
                           'gamma_fast': gf, 'gamma_fast_err': gf_err,
                           'gamma_slow': gs, 'gamma_slow_err': gs_err})
        print(f'  gamma {elapsed(fp):5.0f} s: gamma_fast={gf:.3f}+/-{gf_err:.3f} '
              f'gamma_slow={gs:.3f}+/-{gs_err:.3f}')

# Q ticks: label the mantissa (4..8); the x10^-3 factor is folded into the
# shared x-axis label on the bottom row instead of a separate corner text.
_q_major = [4e-3, 5e-3, 6e-3, 7e-3, 8e-3]
# integer multiples of 10^-3 only; see the note on Figure 3a above
_q_minor = [9e-3]
for a in (axf, axs):
    a.set_xscale('log')
    a.set_yscale('log')
    a.set_xlim(3.25e-3, 9.6e-3)          # margin either side of Q = 3.76-8.27
    add_minor_grid(a)                              # sets minorticks_on + grids
    a.xaxis.set_major_locator(FixedLocator(_q_major))
    a.xaxis.set_minor_locator(FixedLocator(_q_minor))
    a.xaxis.set_major_formatter(FixedFormatter(['4', '5', '6', '7', '8']))
    a.xaxis.set_minor_formatter(NullFormatter())
axf.set_ylabel(r'$\tau_{\mathrm{fast}}$ (s)')
axf.tick_params(labelbottom=False)
axs.set_ylabel(r'$\tau_{\mathrm{slow}}$ (s)')
axs.set_xlabel(r'$Q$ ($\times 10^{-3}\ \AA^{-1}$)')

# --- (d) power-law scaling exponents gamma_fast, gamma_slow vs elapsed time ---
gamma_rows.sort(key=lambda r: r['elapsed'])
xg = [r['elapsed'] for r in gamma_rows]
axg.plot(xg, [r['gamma_fast'] for r in gamma_rows], '-', color='0.75', lw=LW_THIN, zorder=1)
axg.plot(xg, [r['gamma_slow'] for r in gamma_rows], '-', color='0.75', lw=LW_THIN, zorder=1)
for r in gamma_rows:                                 # gamma_fast = circle, gamma_slow = square
    axg.errorbar(r['elapsed'], r['gamma_fast'], yerr=r['gamma_fast_err'], marker='o',
                 color=ecolor(r['fp']), markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
    axg.errorbar(r['elapsed'], r['gamma_slow'], yerr=r['gamma_slow_err'], marker='s',
                 color=ecolor(r['fp']), markersize=MS_SPARSE, capsize=1.5, mfc='none', mew=MEW,
                 elinewidth=LW_THIN, capthick=LW_THIN, zorder=2)
g_handles = [Line2D([], [], marker='o', ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=r'$\gamma_{\mathrm{fast}}$'),
             Line2D([], [], marker='s', ls='none', mfc='none', mec='k', mew=MEW,
                    markersize=MS_SPARSE, label=r'$\gamma_{\mathrm{slow}}$')]
axg.legend(handles=g_handles, loc='lower left')
axg.set_xlabel('Elapsed Time (s)')
axg.set_ylabel(r'Scaling exponent ($\gamma$)')
add_minor_grid(axg)

# (a) and (d) share the elapsed-time axis; widen it once so the first and last
# points are not sitting on the frame (matplotlib's 5 % default is too tight
# for five widely-spaced groups).
_t_pad2 = 0.16 * (max(xe) - min(xe))
axg.set_xlim(min(xe) - _t_pad2, max(xe) + _t_pad2)

# matplotlib's default 5 % y margin leaves the extreme error-bar caps ~1 mm off
# the frame in the autoscaled panels; widen it for the same breathing room the
# explicitly-limited panels have.  (axp sets its own limits above.)
for a in (axf, axs, axg):
    a.set_ymargin(0.11)
    a.autoscale_view()

fig2.tight_layout(pad=0.4, w_pad=1.4, h_pad=0.8)
save_fig(fig2, 'FigureS9_Fit_Parameters.pdf')

# ============================================================
# FIGURE S3 (2 panels): ion-chamber -> photon calibration
#   (left)  air transmission across the calibration measurements (mean ~0.868)
#   (right) photon counts vs upstream IC, with the linear fit CAL_A/CAL_B
# The same CAL_A, CAL_B and AIR_TRANSMISSION feed abs_xsec_coef(), so this
# figure documents exactly the calibration the SAXS panels are scaled by.
# ============================================================
fig3, (axt, axc) = plt.subplots(1, 2, figsize=FIG3_SIZE)
label_panels((axt, axc))

# --- left: air transmission (mean +/- standard deviation reported in legend) ---
axt.plot(np.arange(len(AIR_TRANS_SERIES)), AIR_TRANS_SERIES, 'ko:', ms=MS, mew=MEW, lw=LW_THIN)
axt.axhline(AIR_TRANSMISSION, color='b', ls='-', lw=LW_DATA,
            label=f'mean = {AIR_TRANSMISSION:.4f} ± {AIR_TRANS_STD:.4f}')
axt.set_xlabel('Measurement #')
axt.set_ylabel('Transmission of Air')
# bottom-left: the transmission falls to the right, so that corner is empty
axt.legend(loc='lower left', borderaxespad=0.4)
add_minor_grid(axt)

# --- right: ion chamber -> photon linear calibration ---
_up = CAL_UPIC[CAL_CROP:]
axc.plot(_up, CAL_PHOTONS[CAL_CROP:], 'ko', ms=MS, mew=MEW, label='measured')
_uline = np.linspace(_up.min(), _up.max(), 100)
axc.plot(_uline, CAL_A * _uline + CAL_B, 'b-', lw=LW_DATA,
         label=f'fit: {CAL_A:.2e}·Up_IC {CAL_B:+.2e}')
axc.set_xlabel('Upstream Ion Chamber')
axc.set_ylabel('Number of Photons')
axc.legend(loc='upper left')
add_minor_grid(axc)

fig3.tight_layout(pad=0.4, w_pad=1.4)
save_fig(fig3, 'FigureS3_Calibration.pdf')

plt.show()

"""Figure S6: thermal-cycling repeatability of (VPAVG)30 at 8-ID-I.

This is the direct reversibility control for the SA-XPCS sample: ONE aliquot, in
the SAME aluminium cell and the SAME Quantum Northwest holder used for the
isothermal experiment of Figure 3, taken through seven consecutive 6 -> 34 -> 6 C
cycles.  Unlike the 12-ID-B comparison of Figure 2, which sets two separate
aliquots against each other, every curve here comes from the same material, so a
difference between cycles could only be irreversibility.

Layout (double column):
  top    - the measured sample temperature through the whole 5 h 47 min sequence
  bottom - left:  SAXS I(Q) at the three states sampled in every cycle
           right: g2 at the lowest Q for the two high-temperature states

Temperatures come from /entry/sample/qnw1_temperature.  The QNW stage has three
independently controlled zones holding nine cells; the sample letter selects the
zone (A-C -> qnw1, D-F -> qnw2, G-I -> qnw3), so the B-series sample sits in
qnw1 and the D-series buffer in qnw2.  qnw2 reads 6.00 C throughout, i.e. the
buffer was held cold for the entire sequence.

Acquisition timestamps are NOT usable from the result files: a 2025 reprocessing
overwrote /entry/start_time with the reprocessing date.  They are recovered from
timelist_2022-1.txt, the directory listing of the raw acquisitions, which is the
same source average_ranges.py uses.

Run sequence (each 2 s acquisition at a previously unexposed position):
    B0075 6 C  ->  B0076 ramp to 34 C  ->  D0077 buffer, cooling
    B0078 6 C  ->  B0079 ramp to 34 C  ->  D0080 buffer, cooling
    B0081/83/85/87/89 6 C, each followed by a ramp B0082/84/86/88/90.
The ramp is 270 acquisitions from 6.13 to 34.00 C in 27.7 min = 1.01 C/min; the
cool-down is unmonitored (no acquisitions), so the temperature trace is dashed
across those gaps.

Absolute scale and background subtraction follow Figure 3 exactly: each
acquisition is put on an absolute differential cross section with its OWN
ion-chamber reading via abs_xsec_coef(), the group is averaged, and the averaged
buffer (D0077 + D0080, 89 acquisitions) is subtracted.
"""

import os
import re
import sys
from datetime import datetime

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)
from abs_xsec import abs_xsec_coef, INV_MM_TO_INV_CM
from xpcs_fit import CONTRAST

RESULTS = '/home/8-id-i/2022-1/babnigg202203_nexus/reprocess_results'
TIMELIST = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'timelist_2022-1.txt')

FIG_SIZE = (DOUBLE_COL, 4.75)         # temperature row + the two data panels

# cycle -> (6 C run, ramp run).  The seventh cycle's cool-down was not measured.
CYCLES = [('B0075', 'B0076'), ('B0078', 'B0079'), ('B0081', 'B0082'),
          ('B0083', 'B0084'), ('B0085', 'B0086'), ('B0087', 'B0088'),
          ('B0089', 'B0090')]
BUFFERS = ['D0077', 'D0080']
# the two windows near the top of each ramp, and the mean sample temperature in
# each (measured, not nominal -- see the print-out at the end)
HOT_RANGES = [(241, 250), (261, 270)]
# Colour encodes the STATE (three Okabe-Ito hues, safe under the common colour
# vision deficiencies and ordered cold -> hot); lightness within a hue encodes
# the cycle, light for cycle 1 to dark for cycle 7.  Encoding the cycle as a
# colour in its own right, as a single 7-step map did, made the three states
# impossible to read as groups -- which is the whole point of the panel.
STATE_HUE = ['#0072B2', '#009E73', '#D55E00']    # 6 C, ~32 C, ~34 C
STATE_MARK = ['o', 's', '^']
# The buffer is subtracted with this scale factor.  1.0 leaves the 6 C profile
# with a high-q log-log slope of -0.94; 1.2 brings it to -1.00 with no negative
# points anywhere, i.e. it removes a slight under-subtraction of the solvent.
# (The 12-ID-B reduction uses the same kind of empirical factor, 0.95, on WAXS.)
BG_SCALE = 1.2

SAXS_PATH       = '/xpcs/temporal_mean/scattering_1d'
STATIC_MAP_PATH = '/xpcs/qmap/static_index_mapping'
STATIC_Q_PATH   = '/xpcs/qmap/static_v_list_dim0'
STATIC_PHI_PATH = '/xpcs/qmap/static_v_list_dim1'
TEMP_PATH       = '/entry/sample/qnw1_temperature'
FRAME_TIME_PATH = '/entry/instrument/detector_1/frame_time'
DELAY_PATH      = '/xpcs/multitau/delay_list'
G2_PATH         = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH     = '/xpcs/multitau/normalized_g2_err'
DYN_Q_PATH      = '/xpcs/qmap/dynamic_v_list_dim0'


def load_timelist(path):
    """run-directory name -> acquisition datetime, from the raw `ls -l` dump."""
    out = {}
    for line in open(path):
        p = line.split()
        if len(p) >= 8:
            try:
                out[p[-1]] = datetime.strptime(f'{p[5]} {p[6]}', '%Y-%m-%d %H:%M:%S')
            except ValueError:
                pass
    return out


TIMES = load_timelist(TIMELIST)
_stem = {}                                   # 'B0075' -> full run-name stem
for name in TIMES:
    m = re.match(r'([A-Z]\d{4})_(.*)_Rq0_(\d{5})$', name)
    if m:
        _stem[m.group(1)] = m.group(2)


def run_file(header, n):
    return os.path.join(RESULTS, f'{header}_{_stem[header]}_Rq0_{n:05d}_results.hdf')


def run_time(header, n):
    return TIMES.get(f'{header}_{_stem[header]}_Rq0_{n:05d}')


def n_runs(header):
    return sum(1 for k in TIMES if k.startswith(f'{header}_{_stem[header]}_Rq0_'))


def read_saxs_iq(hf):
    """phi-averaged static I(q), identical to saxpcs.py."""
    inten = np.asarray(hf[SAXS_PATH][()]).reshape(-1)
    idx = hf[STATIC_MAP_PATH][()]
    q_list = hf[STATIC_Q_PATH][()]
    n_phi = hf[STATIC_PHI_PATH].shape[0]
    q_idx = idx // n_phi
    uq = np.unique(q_idx)
    return q_list[uq], np.array([np.nanmean(inten[q_idx == qi]) for qi in uq])


def mean_abs_iq(header, runs):
    """Mean absolute-scale I(Q) [mm^-1] over a list of run numbers.

    Each acquisition is scaled by its own abs_xsec_coef() before averaging, so a
    drift in incident flux across the group is handled per acquisition -- the
    same convention as Figure 3.
    """
    acc, q, n = None, None, 0
    for r in runs:
        f = run_file(header, r)
        if not os.path.exists(f):
            continue
        with h5py.File(f, 'r') as hf:
            q, I = read_saxs_iq(hf)
            I = abs_xsec_coef(hf) * I
        acc = I if acc is None else acc + I
        n += 1
    return q, acc / n, n


def mean_g2(header, runs, q_idx=0):
    """Mean g2 and propagated error over a list of runs, at one dynamic q bin."""
    G, E, tau, qv = [], [], None, None
    for r in runs:
        f = run_file(header, r)
        if not os.path.exists(f):
            continue
        with h5py.File(f, 'r') as hf:
            ft = float(np.asarray(hf[FRAME_TIME_PATH][()]).reshape(-1)[0])
            t = np.asarray(hf[DELAY_PATH][()])
            tau = (t[:, 0] if t.ndim > 1 else t) * ft
            qv = hf[DYN_Q_PATH][()]
            G.append(hf[G2_PATH][()][:, q_idx])
            E.append(hf[G2_ERR_PATH][()][:, q_idx])
    G, E = np.array(G), np.array(E)
    return tau, np.nanmean(G, axis=0), np.sqrt(np.nansum(E**2, axis=0)) / len(G), float(qv[q_idx])


def shades(hue, n):
    """n tints of one hue, light (cycle 1) to dark (cycle n)."""
    base = np.array(mcolors.to_rgb(hue))
    return [tuple((1 - f) * (0.45 + 0.55 * base) + f * (0.55 * base))
            for f in np.linspace(0.0, 1.0, n)]


def read_temp(header, n):
    with h5py.File(run_file(header, n), 'r') as hf:
        return float(np.asarray(hf[TEMP_PATH][()]).reshape(-1)[0])


# ---------------------------------------------------------------- temperature
print('reading temperatures ...')
t0 = run_time('B0075', 1)
trace = []                       # (hours since start, T, header)
for header in [h for c in CYCLES for h in c] + BUFFERS:
    for n in range(1, n_runs(header) + 1):
        tt = run_time(header, n)
        if tt is None or not os.path.exists(run_file(header, n)):
            continue
        trace.append(((tt - t0).total_seconds() / 3600.0, read_temp(header, n), header))
trace.sort()
th = np.array([p[0] for p in trace])
tT = np.array([p[1] for p in trace])
print(f'  {len(trace)} acquisitions, {th[-1]:.2f} h, T = {tT.min():.2f}-{tT.max():.2f} C')

# ------------------------------------------------------------------- buffer
print('averaging buffer ...')
bq, bI, bn = None, None, 0
for b in BUFFERS:
    q, I, n = mean_abs_iq(b, range(1, n_runs(b) + 1))
    bq = q
    bI = I * n if bI is None else bI + I * n
    bn += n
bI /= bn
print(f'  buffer: {bn} acquisitions from {"+".join(BUFFERS)}')

# --------------------------------------------------------- SAXS + XPCS groups
apply_style()
fig = plt.figure(figsize=FIG_SIZE)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.55],
                      hspace=0.42, wspace=0.26,
                      left=0.085, right=0.985, top=0.955, bottom=0.115)
axT = fig.add_subplot(gs[0, :])
axS = fig.add_subplot(gs[1, 0])
axG = fig.add_subplot(gs[1, 1])
label_panels((axT, axS, axG))

SHADES = [shades(h, len(CYCLES)) for h in STATE_HUE]   # [state][cycle]

# temperature: a continuous line inside each measured run, dashed across the
# unmonitored cool-downs (no acquisitions there, so nothing was recorded)
for header in [h for c in CYCLES for h in c] + BUFFERS:
    seg = [(h, T) for h, T, hh in trace if hh == header]
    if seg:
        axT.plot([s[0] for s in seg], [s[1] for s in seg], '-', color='0.35', lw=LW_DATA)


print('averaging SAXS and g2 ...')
saxs_hot, g2_hot = [], []
for ci, (cold, ramp) in enumerate(CYCLES):
    q, I, n = mean_abs_iq(cold, range(1, n_runs(cold) + 1))
    Ia = (I - BG_SCALE * bI) * INV_MM_TO_INV_CM
    pos = Ia > 0
    axS.plot(q[pos], Ia[pos], color=SHADES[0][ci], marker=STATE_MARK[0], ls='none',
             ms=MS_DENSE, mfc='none', mew=LW_THIN)
    tsel = [read_temp(cold, r) for r in (1, n)]
    print(f'  cycle {ci+1} {cold}: {n} acq, T = {np.mean(tsel):.2f} C')

    for hi, (lo, hi_) in enumerate(HOT_RANGES):
        runs = range(lo, hi_ + 1)
        q, I, n = mean_abs_iq(ramp, runs)
        Ia = (I - BG_SCALE * bI) * INV_MM_TO_INV_CM
        pos = Ia > 0
        axS.plot(q[pos], Ia[pos], color=SHADES[hi + 1][ci], marker=STATE_MARK[hi + 1],
                 ls='none', ms=MS_DENSE, mfc='none', mew=LW_THIN)
        Tm = np.mean([read_temp(ramp, r) for r in runs])
        saxs_hot.append(Tm)

        tau, g2, g2e, qv = mean_g2(ramp, runs, q_idx=0)
        m = tau > 0
        axG.errorbar(tau[m], g2[m], yerr=g2e[m], marker=STATE_MARK[hi + 1], ls='none',
                     color=SHADES[hi + 1][ci], ms=MS, mfc='none', mew=LW_THIN,
                     capsize=1.5, elinewidth=LW_THIN, capthick=LW_THIN, alpha=0.9)
        g2_hot.append((Tm, qv))
T_MID = np.mean(saxs_hot[0::2])
T_TOP = np.mean(saxs_hot[1::2])
Q_G2 = g2_hot[0][1]
print(f'  hot windows: runs {HOT_RANGES[0]} -> {T_MID:.1f} C, runs {HOT_RANGES[1]} -> {T_TOP:.1f} C')

# ------------------------------------------------------------------ formatting
axT.set_xlabel('Time (h)')
axT.set_ylabel('Sample temperature (°C)')
axT.set_ylim(0, 38)
add_minor_grid(axT)
# mark the three states that the lower panels sample
for y, hue in ((6, STATE_HUE[0]), (T_MID, STATE_HUE[1]), (T_TOP, STATE_HUE[2])):
    axT.axhline(y, color=hue, ls='--', lw=LW_THIN, zorder=0)
# NOTHING is drawn across the cool-downs: no acquisitions were taken there, so
# the temperature is genuinely unknown and any connecting line would invent it.

axS.set_xscale('log')
axS.set_yscale('log')
axS.set_xlabel(r'$Q$ ($\AA^{-1}$)')
axS.set_ylabel(r'$I(Q)$ (cm$^{-1}$)')
add_minor_grid(axS)
# 3 columns x 2 rows: one column per state, the light swatch for cycle 1 above
# the dark swatch for cycle 7, so the key defines colour AND lightness at once.
labels = ['6 °C', f'{T_MID:.0f} °C', f'{T_TOP:.0f} °C']
handles = []
for si in range(3):
    for ci, tag in ((0, '1'), (len(CYCLES) - 1, str(len(CYCLES)))):
        handles.append(Line2D([], [], marker=STATE_MARK[si], ls='none', mfc='none',
                              mec=SHADES[si][ci], mew=MEW, ms=MS,
                              label=f'{labels[si]}, {tag}'))
leg = axS.legend(handles=handles, loc='upper right', ncol=3, borderaxespad=0.3,
                 title='State, cycle', handlelength=0.8, handletextpad=0.25,
                 columnspacing=0.6, labelspacing=0.2, borderpad=0.3)
leg.get_frame().set_linewidth(LW_THIN)
# headroom so that key sits above the highest curve rather than on it
axS.set_ylim(top=axS.get_ylim()[1] * 9)

axG.set_xscale('log')
# The shortest delays average only 10 acquisitions (against 63-200 in Figure 3),
# so their error bars run far outside the correlation itself; the limits show the
# correlation and let those few caps clip.
axG.set_ylim(0.94, 1.28)
axG.set_xlabel(r'Delay Time, $\tau$ (s)')
axG.set_ylabel('$g_2$')
axG.axhline(1 + CONTRAST, color='0.6', ls=':', lw=LW_THIN)
axG.text(0.985, 1 + CONTRAST, f'$1+\\beta$', transform=axG.get_yaxis_transform(),
         ha='right', va='bottom', color='0.45')
add_minor_grid(axG)
axG.set_title(f'$Q = {Q_G2:.5f}\\ \\AA^{{-1}}$')

save_fig(fig, 'FigureS6_Thermal_Cycle.pdf')
plt.show()

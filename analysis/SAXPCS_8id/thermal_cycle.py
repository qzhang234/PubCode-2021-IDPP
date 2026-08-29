"""Figure S6: thermal-cycling repeatability of (VPAVG)30 at 8-ID-I.

This is the direct reversibility control for the SA-XPCS experiment: ONE aliquot
taken through seven consecutive 6 -> 34 -> 6 C cycles.  Unlike the 12-ID-B
comparison of Figure 2, which sets two separate aliquots against each other,
every curve here comes from the same material, so a difference between cycles
could only be irreversibility.  That is the whole point of the figure and it is
what makes it the more direct test.

It is NOT the same loading as Figure 3.  This is sample S1, measured 2022-03-04;
Figure 3 is sample S3_7, measured 2022-03-06.  Both are aliquots of the same
stock measured in cell position B of the same nine-position Quantum Northwest
holder during the same beamtime, but position B was reloaded in between.  Do not
write "the same aliquot as Figure 3" anywhere.

One consequence is visible in panel (b): the 6 C profile here falls off more
steeply at the lowest q (log-log slope -1.0 below 0.01 A^-1) than the 6 C
reference of Figure 3a (-0.5).  Both are reduced the same way -- same absolute
scaling, same automatic outlier removal (plus the hand cut on cycle 4 described
below), same unscaled buffer subtraction -- and between
0.005 and 0.032 A^-1 the two subtracted profiles agree to within 18 %.  The
difference is in the acquisitions themselves.  Within a single 6 C group here,
the per-acquisition intensity at 0.0035 A^-1 is strongly right-skewed: B0075 has
an RSD of 51 % with a maximum 3.4x its median, B0083 an RSD of 96 % with a
maximum 8.0x its median (skew +2.2 and +3.1).  The SAME acquisitions at 0.03
A^-1 are symmetric with RSD ~ 9 %, i.e. ordinary counting statistics.  That is
what a large object crossing the 10 x 10 um beam in a minority of acquisitions
looks like.  Figure 3's 6 C group does not have it (B0146: RSD 5.8 %, max/median
1.17); the buffers do, D0077 more than D0138.  Part of the difference is a real
property of this loading: it has I(0.004)/I(0.02) = 4.0 against 2.7 for Figure
3's aliquot, and a downstream / upstream ion-chamber ratio of 0.305 against
0.333 (transmissions 0.351 and 0.384 after the air correction), i.e. it was the
less completely dissolved of the two.

MANUAL EXCLUSION -- the one place this figure is not uniformly reduced.
outlier_removal() does NOT catch those low-q spikes.  Its cut is a cosine
similarity on log10 I(q) over the whole q range, so it rejects curves of the
wrong SHAPE; a curve that is 8x high in three low-q bins and normal everywhere
else stays nearly parallel to the group mean and survives.  Eleven acquisitions
of B0083 (cycle 4) -- frames 2, 3, 5, 15, 20, 25, 35, 36, 40, 45 and 47, every
one of them above 1.4x the group median at 0.0035 A^-1 -- are therefore dropped
by hand, through the MANUAL_EXCLUDE table in average_ranges.py, leaving 32 of
50.  That was an author decision taken by inspecting the per-acquisition curves
of the group, not an automatic cut.

  BE CLEAR ABOUT WHAT IT DOES.  Cycle 4 was 31 % ABOVE the mean of the other six
  over 0.004-0.008 A^-1; it is now 11 % BELOW it, and the cycle-to-cycle RSD of
  the 6 C state falls from 12.6 % to 7.7 %.  The cut is applied to cycle 4 ONLY.
  The same 1.4x threshold flags acquisitions in every cycle -- 6, 7, 12, 13, 10,
  7 and 7 for cycles 1-7 -- so applying it uniformly is the fairer comparison.
  Doing that puts cycle 4 at 0.810 against 0.788-0.822 for cycles 1-5, i.e. dead
  centre and no longer an outlier in either direction, and gives an RSD of 9.5 %
  whose residual now comes from cycles 6 and 7 (0.928 and 1.01) rather than from
  cycle 4.  Both numbers are reported in SI Section 7.1.  Switching to the
  uniform variant means extending MANUAL_EXCLUDE to the other six groups and
  re-running ``python average_ranges.py`` for them.

WHY CYCLE 4 STOOD OUT was two different things in the two states, and neither
was irreversibility:

  6 C (circles).  A handful of acquisitions, not the group.  B0083's MEDIAN
  I(0.0035) was 2.80, squarely inside the 2.55-2.96 spread of the other six; its
  MEAN was 4.36 because four of its 43 surviving acquisitions read 22.2, 17.3,
  14.6 and 10.8, i.e. 4-8x the median.  They sat at frames 25, 5, 20 and 45 --
  scattered, not contiguous -- so this was not a drift or a temperature
  excursion but transient objects crossing the beam.  Those four are among the
  eleven now excluded.

  31.9 C (squares).  Not outliers at all: cycle 4's window simply landed hotter.
  The seven ramps cross acquisitions 241-250 at 31.62, 31.87, 31.91, 32.48,
  31.73, 32.28 and 31.68 C, and on the steep part of the transition ln I tracks
  that spread almost perfectly (r = 0.974, d lnI/dT = 0.95 /C).  Removing the
  trend takes the cycle-to-cycle RSD from 32.6 % to 7.1 %.  The 33.8 C window is
  past the steep part and needs no such correction (4.8 % -> 1.4 %).  The
  regression is printed at the end of the run.

Layout (double column):
  top    - the measured sample temperature through the whole 5 h 47 min sequence
  bottom - left:  SAXS I(Q) at the three states sampled in every cycle
           right: g2 at the lowest Q for the two high-temperature states

Temperatures come from /entry/sample/qnw1_temperature.  The QNW stage has three
independently controlled zones holding nine cells; the sample letter selects the
zone (A-C -> qnw1, D-F -> qnw2, G-I -> qnw3), so the B-series sample sits in
qnw1 and the D-series buffer in qnw2.  qnw2 reads 6.00 C throughout, i.e. the
buffer was held cold for the entire sequence.

This script reads nothing from the beamline storage.  Its inputs are the files
average_ranges.py writes into data/: one averaged file per group, carrying the
group's I(q), g2, mean ion-chamber readings and the list of acquisitions that
went into it, and thermal_cycle_temperature.csv, the per-acquisition thermal
history of the whole sequence.

Acquisition timestamps are NOT usable from the result files: a 2025 reprocessing
overwrote /entry/start_time with the reprocessing date.  The elapsed times in
the CSV were recovered by average_ranges.py from timelist_2022-1.txt, the
directory listing of the raw acquisitions.

Run sequence (each 2 s acquisition at a previously unexposed position):
    B0075 6 C  ->  B0076 ramp to 34 C  ->  D0077 buffer, cooling
    B0078 6 C  ->  B0079 ramp to 34 C  ->  D0080 buffer, cooling
    B0081/83/85/87/89 6 C, each followed by a ramp B0082/84/86/88/90.
The ramp is 270 acquisitions from 6.13 to 34.00 C in 27.7 min = 1.01 C/min.
The first two cool-downs WERE recorded, by the buffer runs: D0077 follows the
sample zone from 31.8 to 12.2 C and D0080 from 31.9 to 5.8 C.  No acquisitions
were taken during the other five, so NO TEMPERATURE WAS RECORDED over those
intervals.  They are reconstructed from the control
protocol, which was fixed: cool at 10 C/min until 6 C is reached, then hold at
6 C until the next ramp begins.  Those reconstructed intervals are drawn dashed
and are projections, not measurements; every solid segment is measured.

Absolute scale and background subtraction follow Figure 3 exactly: the group
average is put on an absolute differential cross section by abs_xsec_coef() from
the group's own mean ion-chamber readings, and the averaged
buffer (D0077 + D0080, 78 surviving of 89 acquisitions) is subtracted.  The two
buffer measurements agree to within 6 % at every q below 0.01 A^-1 once filtered
-- I(0.0035) = 1.179 and 1.238 cm^-1, I(0.02) = 0.131 and 0.126 -- and show no
systematic offset above it (only growing bin-to-bin scatter, where both buffers
are weak), so neither is scaled against
the other; scaling only one of them is equivalent to changing BG_SCALE and is
covered by the scan above (D0077 x 2.0 = BG_SCALE 1.44).

EVERY average in this figure runs over the surviving acquisitions and only
those -- the SAXS of each group, its g2, and the mean time and temperature at
which the group is marked in panel (a) -- because all four are read from, or
keyed to, the /xpcs/average/file_list that average_ranges.py wrote alongside the
averaged curves.  The temperature TRACE of panel (a) is the exception, and
deliberately so: it is the thermal history of the experiment, not an average,
and an acquisition whose scattering was rejected still records the temperature
the sample was at.
"""

import glob
import os
import sys
from collections import namedtuple

import numpy as np
import h5py
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from matplotlib.patches import Rectangle
import matplotlib.colors as mcolors

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)), '..'))
from common.acs_style import (DOUBLE_COL, MS, MS_DENSE, MS_SPARSE, MEW, LW_THIN, LW_DATA,
                              apply_style, add_minor_grid, label_panels, save_fig)
from abs_xsec import abs_xsec_coef, INV_MM_TO_INV_CM
from xpcs_fit import CONTRAST

# Everything this figure reads was written by average_ranges.py into data/: one
# averaged file per group, and one CSV holding the per-acquisition thermal
# history of the whole sequence.  Nothing here touches the beamline storage.
DATA = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')
TRACE_CSV = os.path.join(DATA, 'thermal_cycle_temperature.csv')

FIG_SIZE = (DOUBLE_COL, 5.15)         # temperature row, two data panels, key row

# cycle -> (6 C run, ramp run).  The seventh cycle's cool-down was not measured.
CYCLES = [('B0075', 'B0076'), ('B0078', 'B0079'), ('B0081', 'B0082'),
          ('B0083', 'B0084'), ('B0085', 'B0086'), ('B0087', 'B0088'),
          ('B0089', 'B0090')]
BUFFERS = ['D0077', 'D0080']
# Frame ranges, matching the entries in average_ranges.py's FILE_RANGES: the
# whole 6 C group, and two windows near the top of each ramp.  The temperature
# of each window is measured, not nominal -- see the print-out at the end.
COLD_RANGE = (1, 50)
HOT_RANGES = [(241, 250), (261, 270)]
# Colour encodes the CYCLE.  Figure 3's plasma ramp encodes elapsed time, so the
# two figures deliberately do not share a colour mapping.
# Marker SHAPE encodes the temperature state.  Both are keyed along the bottom.
# Seven evenly spaced HUES rather than a sequential ramp: cycle number is a
# label here, not a quantity, and what matters is telling neighbouring curves
# apart.  Raw hsv is unusable for small open markers because its yellow and cyan
# are far lighter than its blue -- they vanish against white -- so each hue is
# darkened to a common perceived luminance (see cycle_colors()).
CYCLE_HUES = np.linspace(0.0, 6.0 / 7.0, 7)
LUMA_CAP = 0.45                      # see cycle_colors()
STATE_MARK = ['o', 's', '^']         # 6 C, ~32 C, ~34 C
# Cooling was measured only twice, by the buffer runs D0077 and D0080 (see the
# docstring).  Nothing was acquired over the other five ramps down; there the
# protocol ran at a fixed rate to 6 C and then held until the next ramp, so the
# gaps can be reconstructed exactly and are drawn dashed as PROJECTIONS.
COOL_RATE = 10.0                     # C/min, fixed
T_HOLD = 6.0                         # C, the temperature held between cycles
# Buffer scale factor.  The 6 C profile still falls at low q rather than showing
# the flat Guinier plateau a fully dispersed solution would, and the obvious
# suspicion is an under-subtracted background.  Scanning the factor shows that
# suspicion is wrong, and the scan is recorded here because the answer is not
# obvious:
#
#   BG    I_6C(0.004)   fraction of the   low-q slope   RSD over the 7 cycles
#                       BG x 1.0 signal   of the 6 C     6 C        34 C
#   1.0      1.411           100 %          -1.04        7.7 %      4.8 %
#   1.2      1.237            88 %          -0.97        8.5 %      4.8 %
#   1.8      0.712            50 %          -0.61       12.2 %      4.8 %
#   2.5      0.101             7 %          +1.03       24.9 %      4.8 %
#
# The buffer is only about 40 % of the 6 C signal at low q (buffer/sample = 0.40
# at 0.0035 A^-1), so flattening the slope means subtracting essentially the
# whole low-q signal: at BG = 2.5 it destroys 93 % of it, turns the slope
# positive, and degrades the cycle-to-cycle repeatability -- the actual result
# of this figure -- from 7.7 % to 25 %, while leaving the 34 C curves untouched
# (4.8 % either way).  A "background" term that changes one curve by 93 % and
# another by 0.4 % is not a shared background.
# The residual low-q rise is therefore treated as sample or cell-window
# scattering, not as mis-subtracted solvent, and is left in.
#
# The factor is 1.0 to match Figures 3 and S8, which subtract the buffer with no
# empirical coefficient at all: there the line is
#     I = coef_sam * I_sample - coef_buf * I_buffer
# with coef_sam and coef_buf the per-acquisition ABSOLUTE-SCALE coefficients
# from abs_xsec_coef(), not a fitted scale.  Using anything but 1.0 here would
# make this figure's background convention differ from the rest of the paper.
BG_SCALE = 1.0

SAXS_PATH       = '/xpcs/temporal_mean/scattering_1d'
STATIC_MAP_PATH = '/xpcs/qmap/static_index_mapping'
STATIC_Q_PATH   = '/xpcs/qmap/static_v_list_dim0'
STATIC_PHI_PATH = '/xpcs/qmap/static_v_list_dim1'
FRAME_TIME_PATH = '/entry/instrument/detector_1/frame_time'
DELAY_PATH      = '/xpcs/multitau/delay_list'
G2_PATH         = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH     = '/xpcs/multitau/normalized_g2_err'
DYN_Q_PATH      = '/xpcs/qmap/dynamic_v_list_dim0'
FILE_LIST       = '/xpcs/average/file_list'


Group = namedtuple('Group', 'q I tau g2 g2e qv hours T n')


def avg_path(header, lo, hi):
    """The averaged file average_ranges.py wrote for one group and frame range."""
    hits = glob.glob(os.path.join(DATA, f'Average_{header}_*_{lo:05d}_{hi:05d}_results.hdf'))
    assert len(hits) == 1, f'expected 1 averaged file for {header} {lo}-{hi}, got {hits}'
    return hits[0]


def load_trace(path):
    """The per-acquisition thermal history written by average_ranges.py.

    Returns the trace -- (hours since the first acquisition, temperature, group
    header) in time order -- and a lookup from dataset name to the same (hours,
    temperature) pair, which is what places each averaged group on panel (a).
    """
    trace, state = [], {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            name, _, secs, T = line.strip().split(',')
            state[name] = (float(secs) / 3600.0, float(T))
            trace.append((*state[name], name.split('_')[0]))
    return trace, state


def read_saxs_iq(hf):
    """phi-averaged static I(q), identical to saxpcs.py."""
    inten = np.asarray(hf[SAXS_PATH][()]).reshape(-1)
    idx = hf[STATIC_MAP_PATH][()]
    q_list = hf[STATIC_Q_PATH][()]
    n_phi = hf[STATIC_PHI_PATH].shape[0]
    q_idx = idx // n_phi
    uq = np.unique(q_idx)
    return q_list[uq], np.array([np.nanmean(inten[q_idx == qi]) for qi in uq])


def group(header, lo, hi):
    """Everything Figure S6 needs from one averaged group.

    q and I(Q) on the absolute scale [mm^-1], the g2 of the lowest q bin with
    its error, the delay times, that bin's q, the mean elapsed time and
    temperature of the group, and how many acquisitions it holds.

    The averaged file was produced by average_ranges.py with the same automatic
    outlier removal, the same cutoff and the same routine applied to the Figure
    3 groups; the ONE exception in this paper is cycle 4 (B0083), which also has
    eleven hand-listed frames removed (see MANUAL EXCLUSION above).  Its
    /xpcs/average/file_list names the acquisitions that survived, and every
    quantity reported here is an average over exactly those: the curves, the g2,
    and the time and temperature at which the group is marked on panel (a).

    The absolute scale enters exactly as it does for Figure 3: the averaged file
    carries the mean incident and transmitted intensities of the same surviving
    acquisitions, and abs_xsec_coef() turns those into one coefficient for the
    averaged curve.  (Until the averaging moved into average_ranges.py this
    figure scaled each acquisition first and then averaged.  The incident flux
    drifts by under 1 % within any of these groups, so the two orders agree to
    0.03 % on average and 0.23 % in the worst q bin of any group -- every number
    this script prints is unchanged to the precision it prints them.)
    """
    with h5py.File(avg_path(header, lo, hi), 'r') as hf:
        q, I = read_saxs_iq(hf)
        I = abs_xsec_coef(hf) * I
        ft = float(np.asarray(hf[FRAME_TIME_PATH][()]).reshape(-1)[0])
        t = np.asarray(hf[DELAY_PATH][()])
        tau = (t[:, 0] if t.ndim > 1 else t) * ft
        qv = float(hf[DYN_Q_PATH][()][0])
        g2 = hf[G2_PATH][()][:, 0]
        g2e = hf[G2_ERR_PATH][()][:, 0]
        kept = [s.decode() if isinstance(s, bytes) else s for s in hf[FILE_LIST][()]]
    kept = [k.replace('_results.hdf', '') for k in kept]
    hrs, T = (float(np.mean([STATE[k][i] for k in kept])) for i in (0, 1))
    return Group(q, I, tau, g2, g2e, qv, hrs, T, len(kept))


def cycle_colors():
    """Seven maximally separated hues, DARKENED to a common luminance ceiling.

    hsv gives the widest hue separation between neighbours but wildly uneven
    lightness: at full value its yellow carries a Rec. 709 luminance of 0.93 and
    its blue only 0.07, so yellow and cyan vanish against white while blue is
    solid.  Only the too-light hues are darkened, by reducing value; the dark
    hues are left alone.  Brightening them instead would clip the RGB channels
    and collapse violet and magenta onto the same colour, which is exactly what
    a symmetric luminance match did.
    """
    out = []
    for h in CYCLE_HUES:
        v = 1.0
        luma = float(np.dot(mcolors.hsv_to_rgb((h, 1.0, 1.0)), (0.2126, 0.7152, 0.0722)))
        if luma > LUMA_CAP:
            v = LUMA_CAP / luma          # darken only; never scale up
        out.append(tuple(mcolors.hsv_to_rgb((h, 1.0, v))))
    return out


# ---------------------------------------------------------------- temperature
print('reading the thermal history ...')
trace, STATE = load_trace(TRACE_CSV)
th = np.array([p[0] for p in trace])
tT = np.array([p[1] for p in trace])
print(f'  {len(trace)} acquisitions, {th[-1]:.2f} h, T = {tT.min():.2f}-{tT.max():.2f} C')

# ------------------------------------------------------------------- buffer
print('averaging buffer ...')
bq, bI, bn = None, None, 0
for b in BUFFERS:
    g = group(b, *COLD_RANGE)
    bq = g.q
    bI = g.I * g.n if bI is None else bI + g.I * g.n
    bn += g.n
bI /= bn
print(f'  buffer: {bn} acquisitions from {"+".join(BUFFERS)}')

# --------------------------------------------------------- SAXS + XPCS groups
apply_style()
fig = plt.figure(figsize=FIG_SIZE)
gs = fig.add_gridspec(2, 2, height_ratios=[1.0, 1.55],
                      hspace=0.42, wspace=0.26,
                      left=0.085, right=0.985, top=0.965, bottom=0.205)
axT = fig.add_subplot(gs[0, :])
axS = fig.add_subplot(gs[1, 0])
axG = fig.add_subplot(gs[1, 1])
label_panels((axT, axS, axG))

COLORS = cycle_colors()

# --- panel (a): measured temperature solid, reconstructed cooling dashed ---
# Each measured run becomes one solid segment, coloured by the cycle it belongs
# to.  Between the end of one measured segment and the start of the next the
# temperature was NOT recorded; it is reconstructed from the known protocol
# (cool at COOL_RATE to T_HOLD, then hold) and drawn dashed.
seg_of = {}
for header in [h for c in CYCLES for h in c] + BUFFERS:
    s = [(h, T) for h, T, hh in trace if hh == header]
    if s:
        seg_of[header] = s
cycle_of = {}
for ci, (cold, ramp) in enumerate(CYCLES):
    cycle_of[cold] = ci
    cycle_of[ramp] = ci
cycle_of['D0077'] = 0
cycle_of['D0080'] = 1

ordered = sorted(seg_of.items(), key=lambda kv: kv[1][0][0])
for header, s in ordered:
    axT.plot([p_[0] for p_ in s], [p_[1] for p_ in s], '-',
             color=COLORS[cycle_of[header]], lw=LW_DATA, zorder=3)
for (h0, s0), (h1, s1) in zip(ordered, ordered[1:]):
    t_end, T_end = s0[-1]
    t_next = s1[0][0]
    if t_next <= t_end:
        continue
    t_cool = t_end + max(T_end - T_HOLD, 0.0) / COOL_RATE / 60.0   # h
    xs = [t_end, min(t_cool, t_next)] + ([t_next] if t_cool < t_next else [])
    ys = [T_end, T_HOLD if t_cool <= t_next else
          T_end - (t_next - t_end) * 60.0 * COOL_RATE] + ([T_HOLD] if t_cool < t_next else [])
    axT.plot(xs, ys, '--', color=COLORS[cycle_of[h0]], lw=LW_THIN, zorder=2)


print('averaging SAXS and g2 ...')
saxs_hot, g2_hot = [], []
# repeatability band: the SI quotes the cycle-to-cycle RSD of the mean
# intensity over 0.004 < q < 0.008 A^-1 for each of the three states
RSD_BAND = (0.004, 0.008)
band_I = [[], [], []]
for ci, (cold, ramp) in enumerate(CYCLES):
    g = group(cold, *COLD_RANGE)
    Ia = (g.I - BG_SCALE * bI) * INV_MM_TO_INV_CM
    pos = Ia > 0
    axS.plot(g.q[pos], Ia[pos], color=COLORS[ci], marker=STATE_MARK[0], ls='none',
             ms=MS_DENSE, mfc='none', mew=LW_THIN)
    band = (g.q >= RSD_BAND[0]) & (g.q <= RSD_BAND[1])
    band_I[0].append(np.nanmean(Ia[band]))
    print(f'  cycle {ci+1} {cold}: {g.n} acq, T = {g.T:.2f} C')
    # mark on (a) exactly where this group was taken
    axT.plot(g.hours, g.T, STATE_MARK[0], color=COLORS[ci],
             mfc='none', mew=MEW, ms=MS_SPARSE, zorder=5)

    for hi, rng in enumerate(HOT_RANGES):
        g = group(ramp, *rng)
        Ia = (g.I - BG_SCALE * bI) * INV_MM_TO_INV_CM
        pos = Ia > 0
        axS.plot(g.q[pos], Ia[pos], color=COLORS[ci], marker=STATE_MARK[hi + 1],
                 ls='none', ms=MS_DENSE, mfc='none', mew=LW_THIN)
        band = (g.q >= RSD_BAND[0]) & (g.q <= RSD_BAND[1])
        band_I[hi + 1].append(np.nanmean(Ia[band]))
        saxs_hot.append(g.T)
        axT.plot(g.hours, g.T, STATE_MARK[hi + 1], color=COLORS[ci],
                 mfc='none', mew=MEW, ms=MS_SPARSE, zorder=5)

        m = g.tau > 0
        axG.errorbar(g.tau[m], g.g2[m], yerr=g.g2e[m], marker=STATE_MARK[hi + 1],
                     ls='none', color=COLORS[ci], ms=MS, mfc='none', mew=LW_THIN,
                     capsize=1.5, elinewidth=LW_THIN, capthick=LW_THIN, alpha=0.9)
        g2_hot.append((g.T, g.qv))
T_MID = np.mean(saxs_hot[0::2])
T_TOP = np.mean(saxs_hot[1::2])
Q_G2 = g2_hot[0][1]
print(f'  hot windows: runs {HOT_RANGES[0]} -> {T_MID:.1f} C, runs {HOT_RANGES[1]} -> {T_TOP:.1f} C')
print(f'  cycle-to-cycle RSD over {RSD_BAND[0]}-{RSD_BAND[1]} A^-1 '
      '(this is what the SI quotes):')
for lab, v in zip((f'{T_HOLD:.0f} C', f'{T_MID:.1f} C', f'{T_TOP:.1f} C'), band_I):
    v = np.array(v)
    print(f'    {lab:>7}: {100 * v.std(ddof=1) / v.mean():5.1f} %   '
          f'per cycle ' + ' '.join(f'{x:.3g}' for x in v))
# The two hot windows are sampled on a ramp, so each cycle crosses them at a
# slightly different temperature.  Near 32 C that is on the steep part of the
# transition and dominates the apparent spread.  Regressing ln I on the MEASURED
# window temperature separates a sampling difference from irreproducibility.
print('  the hot windows are sampled on a ramp, so each cycle crosses them at a '
      'slightly different T:')
for lab, v, Tw in ((f'{T_MID:.1f} C', band_I[1], saxs_hot[0::2]),
                   (f'{T_TOP:.1f} C', band_I[2], saxs_hot[1::2])):
    v, Tw = np.array(v), np.array(Tw)
    sl, ic = np.polyfit(Tw, np.log(v), 1)
    res = np.log(v) - (ic + sl * Tw)
    print(f'    {lab:>7}: windows span {Tw.min():.2f}-{Tw.max():.2f} C, '
          f'dlnI/dT = {sl:.2f} /C, r = {np.corrcoef(Tw, np.log(v))[0, 1]:.3f}  ->  '
          f'RSD {100 * res.std(ddof=1):.1f} % with that trend removed')

# ------------------------------------------------------------------ formatting
axT.set_xlabel('Time (h)')
axT.set_ylabel('Sample temperature (°C)')
axT.set_ylim(0, 38)
add_minor_grid(axT)
# mark the three states that the lower panels sample
for y in (T_HOLD, T_MID, T_TOP):
    axT.axhline(y, color='0.7', ls=':', lw=LW_THIN, zorder=0)

axS.set_xscale('log')
axS.set_yscale('log')
axS.set_xlabel(r'$Q$ ($\AA^{-1}$)')
axS.set_ylabel(r'$I(Q)$ (cm$^{-1}$)')
add_minor_grid(axS)
axG.set_xscale('log')
# The shortest delays average only 10 acquisitions (against 63-200 in Figure 3),
# so their error bars run far outside the correlation itself; the limits show the
# correlation and let those few caps clip.
axG.set_ylim(0.94, 1.28)
axG.set_xlabel(r'Delay Time, $\tau$ (s)')
axG.set_ylabel('$g_2$')
axG.axhline(1 + CONTRAST, color='0.6', ls=':', lw=LW_THIN)
axG.text(0.985, 1 + CONTRAST, '$1+\\beta$', transform=axG.get_yaxis_transform(),
         ha='right', va='bottom', color='0.45')
add_minor_grid(axG)
axG.set_title(f'$Q = {Q_G2:.5f}\\ \\AA^{{-1}}$')

# --- ONE key for the whole figure, drawn by hand ---
# matplotlib's legend can only lay an entry out as handle-then-label, left
# aligned, which cannot produce the layout this key needs: a column header
# CENTRED over its markers, and row labels flush RIGHT against the marker grid.
# The key is therefore drawn directly in figure coordinates -- a 7 x 3 matrix of
# markers, the cycle numbers centred over their columns on the same line as the
# title, the temperature labels right aligned against the grid, and a frame
# around the whole thing.  Everything is sized in POINTS off the one font size,
# so the key scales with the style rather than with the figure.
FS = plt.rcParams['font.size']
ROW_H = 1.40 * FS                    # row pitch
COL_W = 1.65 * FS                    # marker column pitch
GAP   = 1.20 * FS                    # white space between labels and markers
PAD   = 0.55 * FS                    # inside the frame

TITLE = 'Thermal cycle'
ROW_LAB = [f'{T_HOLD:.0f} °C', f'{T_MID:.0f} °C', f'{T_TOP:.0f} °C']

# width of the label column, measured rather than guessed
fig.canvas.draw()
rend = fig.canvas.get_renderer()
probe = [fig.text(0, 0, t) for t in [TITLE] + ROW_LAB]
lab_w = max(t.get_window_extent(rend).width for t in probe) * 72.0 / fig.dpi
for t in probe:
    t.remove()

n_col = len(CYCLES)
box_w = 2 * PAD + lab_w + GAP + n_col * COL_W
box_h = 2 * PAD + 4 * ROW_H
W_pt, H_pt = 72 * fig.get_size_inches()
x0 = (W_pt - box_w) / 2.0            # centred on the figure
y0 = 0.4 * FS                        # just clear of the bottom edge
fx, fy = lambda p: p / W_pt, lambda p: p / H_pt

fig.add_artist(Rectangle((fx(x0), fy(y0)), fx(box_w), fy(box_h),
                         transform=fig.transFigure, facecolor='white',
                         edgecolor='black', lw=LW_THIN, zorder=4))
rows = [y0 + box_h - PAD - (r + 0.5) * ROW_H for r in range(4)]   # top row first
x_lab = x0 + PAD + lab_w                                          # right edge of labels
for r, txt in enumerate([TITLE] + ROW_LAB):
    fig.text(fx(x_lab), fy(rows[r]), txt, ha='right', va='center', zorder=5)
for ci in range(n_col):
    xc = x_lab + GAP + (ci + 0.5) * COL_W
    fig.text(fx(xc), fy(rows[0]), f'{ci + 1}', ha='center', va='center', zorder=5)
    for si in range(3):
        fig.add_artist(Line2D([fx(xc)], [fy(rows[si + 1])], marker=STATE_MARK[si],
                              ls='none', mfc='none', mec=COLORS[ci], mew=MEW,
                              ms=MS, transform=fig.transFigure, zorder=5))

save_fig(fig, 'FigureS6_Thermal_Cycle.pdf')
plt.show()

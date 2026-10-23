"""Reduce the 2021-3 temperature series the same way average_ranges.py reduces 2022-1.

BEAMLINE-ONLY.  Reads /gdata, which mounts on the 10.x hosts (amber) and not on
the 164.x analysis machines.  Everything downstream reads only what this writes
into data_2021-3/.

THE DATA
--------
Four ELP loadings (D6-D9) in the same cell block, each taken from 20 C to a
different final temperature on the QNW zone-3 controller, on 17 Oct 2021:

    H0402  D8  28.00 C        H0407  D6  28.25 C
    H0404  D7  28.50 C        H0413  D9  28.75 C

H0412 (D9 at 28.25 C) was stopped at 21.1 C after 11 minutes and never reached
temperature, so D9's run is H0413.  Every run starts at 20 C and climbs at
0.10 C/min, so the target is reached about 80 minutes in.  That is the main
difference from the 2022-1 series behind Figure 3, which was brought from 20 to
30 C at 1 C/min.

WHAT IS WRITTEN
---------------
1. TRACE: one row per acquisition with its time, zone-3 temperature, the two ion
   chambers and the mean I(Q) over SPIKE_BAND, 0.004-0.008 A^-1, the band
   Figure S4b reports.  No averaging, so the onset can be placed in temperature.
2. AVERAGES: 100-acquisition blocks, run through average_ranges.process_group so
   the two outlier cuts are exactly the ones every 2022-1 group gets.  Blocks are
   aligned so that one boundary falls on the first acquisition at the target
   temperature, which is t_w = 0, the analogue of the RTD reaching 30 C.
3. The 6 C reference of each loading and the 12IDB static standard B0005, each
   averaged whole, for the SAXS baseline and for the speckle contrast.
4. The empty-cell and air measurements, each header averaged whole through the
   same two outlier cuts, as candidate backgrounds.  They were correlated with
   the same qmap and boost_corr flags as the sample data.

    python reduce_babnigg202110.py            # everything
    python reduce_babnigg202110.py empties    # only item 4
"""

import os
import sys
import glob

import numpy as np
import h5py

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import average_ranges as ar          # noqa: E402  the 2022-1 reduction, reused as is

PREFIX = '/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results'
OUT_DIR = os.path.join(HERE, 'data_2021-3')
TRACE_CSV = os.path.join(OUT_DIR, 'trace_2021-3.csv')

RUNS = ['H0402_D8_280C10p', 'H0404_D7_285C10p', 'H0407_D6_282C10p',
        'H0413_D9_287C10p']
REFERENCES = ['H0402_D8_060C10p', 'H0403_D8_060C10p', 'H0404_D7_060C10p',
              'H0407_D6_060C10p', 'H0413_D9_060C10p']
STANDARD = 'B0005_12IDB_Standard_006C'
# Candidate backgrounds.  E0008 and E0304 (9-10 Oct) carry no transmitted-beam
# monitor and an incident monitor on a different scale from 17 Oct, so only the
# 17 Oct ones can go through abs_xsec_coef as recorded.
EMPTIES = ['E0008__006C', 'E0304_Empty_006C', 'E0408_Air_006C',
           'E0411_Empty_3_006C', 'G0409_Empty_006C', 'G0410_Empty_2_006C']
BLOCK = 100                           # acquisitions per average, as in 2022-1 (100-150)
AT_TARGET = 0.01                      # C; zone 3 counts as at temperature within this

# average_ranges.py reads these two module globals, so pointing them here sends
# its unmodified code path at the 2021-3 tree.  Nothing about 2022-1 changes.
ar.prefix = PREFIX
ar.out_dir = OUT_DIR

TEMP = '/entry/sample/qnw3_temperature'
TEMP_SET = '/entry/sample/qnw3_temperature_set'


def scalar(hf, path):
    return float(np.asarray(hf[path][()]).ravel()[0])


def trace(run):
    """Per-acquisition rows for one run: frame, time, T, T_set, I0, I_T, band mean."""
    files = ar.group_files(run)
    q = ar.static_q(files[0])
    band = (q >= ar.SPIKE_BAND[0]) & (q < ar.SPIKE_BAND[1])
    rows = []
    for f in files:
        with h5py.File(f, 'r') as hf:
            saxs = np.asarray(hf[ar.SAXS_PATH][()]).reshape(-1)
            rows.append([ar.extract_frame(f), ar.get_start_time(f),
                         scalar(hf, TEMP), scalar(hf, TEMP_SET),
                         scalar(hf, ar.INCIDENT_PATH), scalar(hf, ar.TRANSMITTED_PATH),
                         float(np.nanmean(saxs[band]))])
    return rows


def blocks(rows):
    """100-acquisition frame ranges with one boundary at the first frame on target."""
    frames = np.array([r[0] for r in rows])
    T, Tset = np.array([r[2] for r in rows]), np.array([r[3] for r in rows])
    a0 = int(frames[np.argmax(T >= Tset - AT_TARGET)])
    first, last = int(frames.min()), int(frames.max())
    # Block starts every BLOCK frames, counted both ways from a0.  The earliest
    # block is whatever is left before the first full one, so it may be short.
    starts = sorted(set([first] + [s for s in range(a0, first - 1, -BLOCK) if s > first]
                        + list(range(a0, last + 1, BLOCK))))
    ends = [s - 1 for s in starts[1:]] + [last]
    return a0, list(zip(starts, ends))


def reduce_whole(group):
    """Average every acquisition of a group into one file."""
    frames = [ar.extract_frame(f) for f in ar.group_files(group)]
    ar.process_group(group, [(min(frames), max(frames))])


def main(only=None):
    os.makedirs(OUT_DIR, exist_ok=True)
    if only == 'empties':
        for group in EMPTIES:
            reduce_whole(group)
        return

    all_rows, ranges = [], {}
    for run in RUNS:
        rows = trace(run)
        a0, rng = blocks(rows)
        ranges[run] = rng
        all_rows += [[run] + r for r in rows]
        print(f'{run}: {len(rows)} acquisitions, on target from frame {a0}, '
              f'{len(rng)} blocks')
    with open(TRACE_CSV, 'w') as fh:
        fh.write('run,frame,start_time,qnw3_C,qnw3_set_C,incident,transmitted,'
                 'I_band\n')
        for r in all_rows:
            fh.write(','.join([r[0], str(r[1]), r[2], f'{r[3]:.3f}', f'{r[4]:.3f}',
                               f'{r[5]:.6g}', f'{r[6]:.6g}', f'{r[7]:.6g}']) + '\n')
    print(f'trace -> {os.path.basename(TRACE_CSV)} ({len(all_rows)} rows)')

    for run in RUNS:
        ar.process_group(run, ranges[run])
    for ref in REFERENCES:
        ar.process_group(ref, [(1, 20)])
    ar.process_group(STANDARD, [(1, 100)])


if __name__ == '__main__':
    main(sys.argv[1] if len(sys.argv) > 1 else None)

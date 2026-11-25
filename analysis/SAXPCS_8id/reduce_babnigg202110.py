"""Write data_2021-3/: the per-acquisition trace and the group averages.

BEAMLINE-ONLY: reads /gdata (amber, adamite, amazonite), and averages through
average_ranges.py, which is not part of this repository.  Everything
downstream reads only what this writes.

    python reduce_babnigg202110.py            # the trace, the groups, the 6 C references, the standard
    python reduce_babnigg202110.py empties    # the six candidate backgrounds
    python reduce_babnigg202110.py ramp20     # the ~20 C reference: first 20 acquisitions of each ramp

Run both after a reprocessing: the first call does not rebuild the background
averages.  The repository README describes what is written.
"""

import os
import sys

import h5py
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, HERE)
import average_ranges as ar          # noqa: E402  the reduction behind every group average

PREFIX = '/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results'
OUT_DIR = os.path.join(HERE, 'data_2021-3')
TRACE_CSV = os.path.join(OUT_DIR, 'trace_2021-3.csv')

RUNS = ['H0402_D8_280C10p', 'H0404_D7_285C10p', 'H0407_D6_282C10p',
        'H0413_D9_287C10p']
REFERENCES = ['H0402_D8_060C10p', 'H0403_D8_060C10p', 'H0404_D7_060C10p',
              'H0407_D6_060C10p', 'H0413_D9_060C10p']
STANDARD = 'B0005_12IDB_Standard_006C'
# Candidate backgrounds.  E0008 and E0304 (9-10 Oct) carry no transmitted-beam
# monitor and an incident monitor on a different scale from 17 Oct.
EMPTIES = ['E0008__006C', 'E0304_Empty_006C', 'E0408_Air_006C',
           'E0411_Empty_3_006C', 'G0409_Empty_006C', 'G0410_Empty_2_006C']
BLOCK = 100                           # acquisitions per group
AT_TARGET = 0.01                      # C; zone 3 counts as at its target within this
# The ~20 C reference of each loading: its first RAMP_REF acquisitions (RTD
# 20.0-20.2 C), the same number for all four.  They go to their own folder so
# the group lists, which glob data_2021-3/, never pick them up.
RAMP_REF = 20
RAMP_DIR = os.path.join(OUT_DIR, 'ramp20')

# average_ranges.py reads these two module variables; pointing them here sends
# its code at the 2021-3 tree.
ar.prefix = PREFIX
ar.out_dir = OUT_DIR

TEMP = '/entry/sample/qnw3_temperature'
TEMP_SET = '/entry/sample/qnw3_temperature_set'


def scalar(hf, path):
    return float(np.asarray(hf[path][()]).ravel()[0])


def trace(run):
    """Per-acquisition rows of one run: frame, start time, T, T_set, I0, I_T,
    and the mean of the stored I(Q) entries over 0.004-0.008 A^-1."""
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
    """Frame ranges of BLOCK acquisitions, with one boundary at the first frame
    on target (t_w = 0).  The earliest range takes what is left before the first
    full one, so it may be short."""
    frames = np.array([r[0] for r in rows])
    T = np.array([r[2] for r in rows])
    T_set = np.array([r[3] for r in rows])
    on_target = T >= T_set - AT_TARGET
    if not on_target.any():
        raise ValueError('the run never reaches its target temperature')
    a0 = int(frames[np.argmax(on_target)])
    first = int(frames.min())
    last = int(frames.max())

    starts = [first]
    s = a0
    while s > first:                  # every BLOCK frames below a0
        starts.append(s)
        s = s - BLOCK
    s = a0
    while s <= last:                  # a0 and every BLOCK frames above it
        starts.append(s)
        s = s + BLOCK
    starts = sorted(set(starts))

    ranges = []
    for k in range(len(starts)):
        if k + 1 < len(starts):
            ranges.append((starts[k], starts[k + 1] - 1))
        else:
            ranges.append((starts[k], last))
    return a0, ranges


def reduce_whole(group):
    """Average every acquisition of a group into one file."""
    frames = [ar.extract_frame(f) for f in ar.group_files(group)]
    ar.process_group(group, [(min(frames), max(frames))])


def main(only):
    os.makedirs(OUT_DIR, exist_ok=True)
    if only == 'empties':
        for group in EMPTIES:
            reduce_whole(group)
        return
    if only == 'ramp20':
        os.makedirs(RAMP_DIR, exist_ok=True)
        ar.out_dir = RAMP_DIR
        for run in RUNS:
            ar.process_group(run, [(1, RAMP_REF)])
        return

    all_rows = []
    ranges = {}
    for run in RUNS:
        rows = trace(run)
        a0, run_ranges = blocks(rows)
        ranges[run] = run_ranges
        for r in rows:
            all_rows.append([run] + r)
        print(f'{run}: {len(rows)} acquisitions, on target from frame {a0}, '
              f'{len(run_ranges)} blocks')
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
    if len(sys.argv) > 1:
        main(sys.argv[1])
    else:
        main(None)

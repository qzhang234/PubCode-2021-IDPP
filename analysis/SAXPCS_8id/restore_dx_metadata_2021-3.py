"""Carry the ion chambers and ring current from the 2021-3 archive into its NeXus tree.

BEAMLINE-ONLY.  Needs /gdata, which mounts on the 10.x hosts (amber) and not on
the 164.x analysis machines.  Kept here for provenance.

WHAT IT DOES
------------
The absolute scattering cross section of abs_xsec.py needs, for every
acquisition, the incident and transmitted beam intensities.  In the 2021-3 Data
Exchange archive (babnigg202110) they are recorded per acquisition as

    NeXus destination                                         Data Exchange source
    /entry/instrument/incident_beam/incident_beam_intensity      source_end/I0Monitor
    /entry/instrument/incident_beam/transmitted_beam_intensity   source_end/TransmissionMonitor
    /entry/instrument/incident_beam/ring_current                 source_end/current

the same mapping restore_dx_metadata.py uses for 2022-1.  This script checks both
levels of the NeXus tree, the per-acquisition *_metadata.hdf and every
reprocess_results/*_results.hdf, against the archive value of the same
acquisition.  It writes only where the NeXus value is missing, not finite, or
differs from the archive, so it is idempotent and safe to rerun.

nexus_convert_babnigg202110.py already wrote these fields when it built the
tree, and boost_corr copies them into each result, so on a tree built by that
script the expected outcome is that nothing needs writing.  This script is the
check that the outcome is what was intended.

WHAT IT CANNOT DO
-----------------
The 9-10 Oct measurements (the 12IDB standard, the 150 nm silica, E0008 and
E0304) recorded no TransmissionMonitor, and their I0Monitor reads near 1e9
rather than the near 9e5 of 17 Oct, so the ion-chamber electronics were set up
differently.  A missing archive field is reported and left as it is; nothing is
invented.

Usage
-----
    python restore_dx_metadata_2021-3.py            # dry run; writes the manifest only
    python restore_dx_metadata_2021-3.py --apply    # perform the writes
"""

import os
import csv
import glob
import argparse
import multiprocessing

import numpy as np
import h5py

ARCHIVE = '/gdata/s8id-dmdtn/2021-3/babnigg202110'
NEXUS = '/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus'
RESULTS = os.path.join(NEXUS, 'reprocess_results')
MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        'dx_restore_manifest_2021-3.csv')

FIELDS = (
    ('/entry/instrument/incident_beam/incident_beam_intensity',
     'measurement/instrument/source_end/I0Monitor'),
    ('/entry/instrument/incident_beam/transmitted_beam_intensity',
     'measurement/instrument/source_end/TransmissionMonitor'),
    ('/entry/instrument/incident_beam/ring_current',
     'measurement/instrument/source_end/current'),
)


def scalar(obj):
    """First element of a field stored as (), (1,) or (1, 1)."""
    return float(np.asarray(obj[()]).ravel()[0])


def archive_index():
    """Acquisition name -> its Data Exchange file, for every acquisition in the archive.

    Most scans are a directory of acquisition sub-directories; a few early ones
    are a single acquisition sitting directly in the scan directory.
    """
    index = {}
    for scan_dir in glob.glob(os.path.join(ARCHIVE, '*')):
        if not os.path.isdir(scan_dir):
            continue
        subs = [d for d in glob.glob(os.path.join(scan_dir, os.path.basename(scan_dir) + '_*'))
                if os.path.isdir(d)]
        for d in subs or [scan_dir]:
            found = sorted(glob.glob(os.path.join(d, '*_0001-*.hdf')))
            if found:
                index[os.path.basename(d)] = found[0]
    return index


def read_archive(path):
    """{nexus_path: value} for the fields this acquisition recorded."""
    out = {}
    try:
        with h5py.File(path, 'r') as hf:
            for nx_path, dx_path in FIELDS:
                if dx_path in hf:
                    out[nx_path] = scalar(hf[dx_path])
    except OSError:
        pass
    return out


def check_file(path, values, apply_changes):
    """Compare one NeXus file with the archive values; write if asked.  Returns rows."""
    rows = []
    with h5py.File(path, 'r+' if apply_changes else 'r') as hf:
        for nx_path, _dx in FIELDS:
            want = values.get(nx_path)
            have = scalar(hf[nx_path]) if nx_path in hf else None
            if want is None:
                rows.append((nx_path, have, None, 'not recorded in archive'))
                continue
            if have is not None and np.isfinite(have) and np.isclose(have, want, rtol=1e-6, atol=0):
                rows.append((nx_path, have, want, 'already correct'))
                continue
            action = 'written' if apply_changes else 'would write'
            if apply_changes:
                if nx_path in hf:
                    ds = hf[nx_path]
                    ds[...] = np.asarray(want, dtype=ds.dtype).reshape(ds.shape)
                else:
                    hf.create_dataset(nx_path, data=float(want))
            rows.append((nx_path, have, want, action))
    return rows


def process(job):
    acq, dx_file, apply_changes = job
    values = read_archive(dx_file)
    targets = [os.path.join(NEXUS, acq, f'{acq}_metadata.hdf'),
               os.path.join(RESULTS, f'{acq}_results.hdf')]
    rows = []
    for t in targets:
        if not os.path.exists(t):
            continue
        try:
            for nx_path, have, want, action in check_file(t, values, apply_changes):
                rows.append([acq, os.path.basename(t), nx_path.split('/')[-1], have, want, action])
        except OSError as exc:
            rows.append([acq, os.path.basename(t), '', None, None, f'unreadable: {exc}'])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--apply', action='store_true', help='perform the writes')
    ap.add_argument('--jobs', type=int, default=24)
    args = ap.parse_args()

    index = archive_index()
    acqs = sorted(d for d in os.listdir(NEXUS)
                  if os.path.isdir(os.path.join(NEXUS, d)) and d in index)
    print(f'{len(acqs)} NeXus acquisitions with an archive record')
    with multiprocessing.Pool(args.jobs) as pool:
        rows = [r for part in pool.map(process, [(a, index[a], args.apply) for a in acqs],
                                       chunksize=16) for r in part]

    # The manifest records every field that was not already correct: what was
    # written, and what the archive never recorded.  The totals above it cover
    # the rest.
    with open(MANIFEST, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['acquisition', 'file', 'field', 'nexus_before', 'archive', 'action'])
        w.writerows(r for r in rows if r[5] != 'already correct')

    files = {(r[0], r[1]) for r in rows}
    print(f'{len(files)} files checked '
          f'({sum(1 for f in files if f[1].endswith("_metadata.hdf"))} metadata, '
          f'{sum(1 for f in files if f[1].endswith("_results.hdf"))} results)')
    tally = {}
    for r in rows:
        tally[(r[2], r[5])] = tally.get((r[2], r[5]), 0) + 1
    for (field, action), n in sorted(tally.items()):
        print(f'  {field:28s} {action:26s} {n:6d}')
    print(f'manifest -> {MANIFEST}')
    if not args.apply and any(r[5] == 'would write' for r in rows):
        print('rerun with --apply to write')


if __name__ == '__main__':
    main()

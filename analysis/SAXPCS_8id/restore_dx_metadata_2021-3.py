"""Check that every file of the 2021-3 NeXus tree carries its ion chambers and
ring current, and copy them from the Data Exchange archive where it does not.

BEAMLINE-ONLY: needs /gdata (amber, adamite, amazonite).

For every acquisition, the per-acquisition *_metadata.hdf and the
reprocess_results/*_results.hdf are compared with the archive value of

    NeXus destination                                         Data Exchange source
    /entry/instrument/incident_beam/incident_beam_intensity      source_end/I0Monitor
    /entry/instrument/incident_beam/transmitted_beam_intensity   source_end/TransmissionMonitor
    /entry/instrument/incident_beam/ring_current                 source_end/current

and written only where the NeXus value is missing, not finite, or different,
so a rerun is safe.  A field the archive never recorded (the 9-10 Oct
measurements have no TransmissionMonitor) is reported and left alone.  Every
field that was not already correct is listed in dx_restore_manifest_2021-3.csv.

    python restore_dx_metadata_2021-3.py            # dry run; writes the manifest only
    python restore_dx_metadata_2021-3.py --apply    # perform the writes
"""

import argparse
import csv
import glob
import multiprocessing
import os

import h5py
import numpy as np

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
        acquisition_dirs = []
        for d in glob.glob(os.path.join(scan_dir, os.path.basename(scan_dir) + '_*')):
            if os.path.isdir(d):
                acquisition_dirs.append(d)
        if len(acquisition_dirs) == 0:
            acquisition_dirs = [scan_dir]
        for d in acquisition_dirs:
            found = sorted(glob.glob(os.path.join(d, '*_0001-*.hdf')))
            if found:
                index[os.path.basename(d)] = found[0]
    return index


def read_archive(path):
    """{nexus_path: value} for the fields this acquisition recorded."""
    values = {}
    with h5py.File(path, 'r') as hf:
        for nx_path, dx_path in FIELDS:
            if dx_path in hf:
                values[nx_path] = scalar(hf[dx_path])
    return values


def check_file(path, values, apply_changes):
    """Compare one NeXus file with the archive values; write if asked.  Returns rows."""
    rows = []
    if apply_changes:
        mode = 'r+'
    else:
        mode = 'r'
    with h5py.File(path, mode) as hf:
        for nx_path, dx_path in FIELDS:
            want = values.get(nx_path)
            if nx_path in hf:
                have = scalar(hf[nx_path])
            else:
                have = None
            if want is None:
                rows.append((nx_path, have, None, 'not recorded in archive'))
                continue
            if have is not None and np.isfinite(have) and np.isclose(have, want, rtol=1e-6, atol=0):
                rows.append((nx_path, have, want, 'already correct'))
                continue
            if apply_changes:
                if nx_path in hf:
                    ds = hf[nx_path]
                    ds[...] = np.asarray(want, dtype=ds.dtype).reshape(ds.shape)
                else:
                    hf.create_dataset(nx_path, data=float(want))
                rows.append((nx_path, have, want, 'written'))
            else:
                rows.append((nx_path, have, want, 'would write'))
    return rows


def process(job):
    """Check both NeXus files of one acquisition.  job = (acquisition, archive file, apply)."""
    acq, dx_file, apply_changes = job
    values = read_archive(dx_file)
    targets = [os.path.join(NEXUS, acq, f'{acq}_metadata.hdf'),
               os.path.join(RESULTS, f'{acq}_results.hdf')]
    rows = []
    for target in targets:
        if not os.path.exists(target):
            continue
        for nx_path, have, want, action in check_file(target, values, apply_changes):
            rows.append([acq, os.path.basename(target), nx_path.split('/')[-1], have, want, action])
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--apply', action='store_true', help='perform the writes')
    ap.add_argument('--jobs', type=int, default=24)
    args = ap.parse_args()

    index = archive_index()
    acqs = []
    for d in os.listdir(NEXUS):
        if os.path.isdir(os.path.join(NEXUS, d)) and d in index:
            acqs.append(d)
    acqs = sorted(acqs)
    print(f'{len(acqs)} NeXus acquisitions with an archive record')
    jobs = []
    for a in acqs:
        jobs.append((a, index[a], args.apply))
    with multiprocessing.Pool(args.jobs) as pool:
        results = pool.map(process, jobs, chunksize=16)
    rows = []
    for part in results:
        rows.extend(part)

    # the manifest lists every field that was not already correct
    with open(MANIFEST, 'w', newline='') as fh:
        w = csv.writer(fh)
        w.writerow(['acquisition', 'file', 'field', 'nexus_before', 'archive', 'action'])
        for r in rows:
            if r[5] != 'already correct':
                w.writerow(r)

    files = set()
    for r in rows:
        files.add((r[0], r[1]))
    n_metadata = 0
    n_results = 0
    for f in files:
        if f[1].endswith('_metadata.hdf'):
            n_metadata += 1
        if f[1].endswith('_results.hdf'):
            n_results += 1
    print(f'{len(files)} files checked ({n_metadata} metadata, {n_results} results)')
    tally = {}
    for r in rows:
        tally[(r[2], r[5])] = tally.get((r[2], r[5]), 0) + 1
    for (field, action), n in sorted(tally.items()):
        print(f'  {field:28s} {action:26s} {n:6d}')
    print(f'manifest -> {MANIFEST}')
    would_write = False
    for r in rows:
        if r[5] == 'would write':
            would_write = True
    if not args.apply and would_write:
        print('rerun with --apply to write')


if __name__ == '__main__':
    main()

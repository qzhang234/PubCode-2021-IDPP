"""Average XPCS results over explicit frame-number ranges.

The groups (file-name headers) and their frame ranges are set inline in
``FILE_RANGES`` below.

For each range this script:
  * globs the group's result files and selects those whose frame number is in
    the range,
  * removes outliers (cross-correlation threshold on g2, g2_err, saxs_1d),
  * averages g2, g2_err and saxs_1d over the surviving files,
  * ALSO averages the two scalar ion-chamber monitors -- the incident
    (upstream) and transmitted (downstream) beam intensities -- over the same
    surviving files.  These are needed downstream in saxpcs.py to compute the
    absolute scattering cross-section coefficient of each averaged dataset, so
    each frame range carries its own mean flux/transmission.  The scalars are
    averaged SEPARATELY from g2/saxs_1d and are deliberately NOT put through
    outlier_removal (that routine cross-correlates every field it is handed and
    is meaningless for single scalars),
  * writes the averaged data into a new HDF file that keeps the same structure
    as the originals (a copy of the first file in the range is used as the
    template), inside an ``average`` folder one level above ``reprocess_results``.

The averaged file name starts with 'Average_' and carries the frame range so it
is easy to tell apart from the raw files, e.g.::

    Average_B0147_S3_7_300C10p_att00_Rq0_00950_01050_results.hdf

In the averaged file:
  * ``/entry/start_time`` is set to the acquisition time of the first file in
    the range (looked up in the beamline time list),
  * ``/xpcs/average/file_list`` is added, listing every file included in the
    average (i.e. the files that survived outlier removal).

The plotting is done separately in ``Plot_Ave_Ranges.py``.
"""

import sys
import os
import re
import glob
import shutil

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common.utils import outlier_removal, average_datasets

# --- PARAMETERS ---
prefix = '/home/8-id-i/2022-1/babnigg202203_nexus/reprocess_results'

# Groups to average: file-name header -> list of (start, end) frame ranges
# (inclusive), one averaged file per range.  Usually one range per group.

FILE_RANGES = {
    'B0147': [(1, 200),
        (200, 350),
        (351, 500),
        (501, 650),
        (651, 800),
        (801, 950),
        (951, 1050),
        (1051, 1150),
        (1151, 1250),
        (1251, 1313)],
    'B0146': [(1, 50)],
    'D0138': [(1, 50)],
}


# Averaged files are written into the local 'data' folder next to this script
# (created if it does not already exist), so downstream analysis reads locally
# instead of from the remote /home/8-id-i beamline storage.
out_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), 'data')

# Time list used to recover the acquisition time of each raw dataset.
timelist_path = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                             'timelist_2022-1.txt')

# Percentile cutoff for the cross-correlation outlier removal.
percentile = 5

# --- HDF field locations (same layout in the raw and averaged files) ---
G2_PATH     = '/xpcs/multitau/normalized_g2'
G2_ERR_PATH = '/xpcs/multitau/normalized_g2_err'
SAXS_PATH   = '/xpcs/temporal_mean/scattering_1d'
START_TIME  = '/entry/start_time'
FILE_LIST   = '/xpcs/average/file_list'
# Scalar ion-chamber monitors, populated from the raw cluster results.  These
# are averaged per range (separately from the array fields) and written back so
# saxpcs.py can derive an absolute cross-section coefficient per averaged file.
INCIDENT_PATH    = '/entry/instrument/incident_beam/incident_beam_intensity'
TRANSMITTED_PATH = '/entry/instrument/incident_beam/transmitted_beam_intensity'

_frame_re = re.compile(r'_(\d+)_results')


def extract_frame(fname):
    """Return the integer frame number embedded in a result file name."""
    m = _frame_re.search(os.path.basename(fname))
    return int(m.group(1)) if m else -1


def load_timelist(path):
    """Parse an ``ls -l`` style time list into ``{folder_name: 'YYYY-MM-DD HH:MM:SS'}``.

    The dataset folder name matches a result file name with ``_results.hdf``
    stripped off.
    """
    times = {}
    with open(path) as fh:
        for line in fh:
            tokens = line.split()
            # dirs/files: perms links owner group size date time name
            if len(tokens) < 8:
                continue
            name = tokens[-1]
            date = tokens[-3]
            clock = tokens[-2]
            times[name] = f'{date} {clock}'
    return times


def get_start_time(fname, timelist):
    """Look up the acquisition time for a result file, or None if absent."""
    key = os.path.basename(fname).replace('_results.hdf', '')
    return timelist.get(key)


def read_range_data(section_files):
    """Read g2, g2_err, saxs_1d and the two beam monitors for every file.

    Returns
    -------
    data_dict : dict of stacked arrays (first axis = file index) for the fields
        that go through outlier removal (g2, g2_err, saxs_1d).
    read_files : list of the files read successfully, in stack order.
    incident : 1d array of the incident (upstream) beam intensity per file.
    transmitted : 1d array of the transmitted (downstream) beam intensity per file.

    ``incident``/``transmitted`` are kept OUT of ``data_dict`` on purpose: they
    are scalars per file and must not be fed to ``outlier_removal`` (which
    cross-correlates every key it receives).  They stay aligned with
    ``read_files`` so the g2-based outlier mask can be applied to them too.
    """
    data_dict = {'g2': [], 'g2_err': [], 'saxs_1d': []}
    incident, transmitted = [], []
    read_files = []
    for f in section_files:
        try:
            with h5py.File(f, 'r') as hf:
                g2 = hf[G2_PATH][()]
                g2_err = hf[G2_ERR_PATH][()]
                saxs = np.asarray(hf[SAXS_PATH][()])
                # collapse a (1, N)/segmented SAXS to a single 1d curve
                if saxs.ndim == 2:
                    saxs = np.nanmean(saxs, axis=0)
                inc = float(np.asarray(hf[INCIDENT_PATH][()]).ravel()[0])
                trn = float(np.asarray(hf[TRANSMITTED_PATH][()]).ravel()[0])
        except Exception as exc:
            print(f'  failed to read {os.path.basename(f)}: {exc}')
            continue
        data_dict['g2'].append(g2)
        data_dict['g2_err'].append(g2_err)
        data_dict['saxs_1d'].append(saxs)
        incident.append(inc)
        transmitted.append(trn)
        read_files.append(f)

    for k in data_dict:
        data_dict[k] = np.array(data_dict[k])
    return data_dict, read_files, np.array(incident), np.array(transmitted)


def _write_scalar(hf, path, value):
    """Overwrite an existing scalar dataset, keeping its stored shape/dtype."""
    ds = hf[path]
    ds[...] = np.asarray(value, dtype=ds.dtype).reshape(ds.shape)


def save_average(template, out_path, avg_dict, start_time, included_files,
                 incident_avg, transmitted_avg):
    """Copy ``template`` to ``out_path`` and overwrite it with the averaged data."""
    shutil.copyfile(template, out_path)

    with h5py.File(out_path, 'r+') as hf:
        hf[G2_PATH][...]     = avg_dict['g2']
        hf[G2_ERR_PATH][...] = avg_dict['g2_err']
        # write back with the file's own SAXS shape (e.g. (1, N))
        hf[SAXS_PATH][...] = avg_dict['saxs_1d'].reshape(hf[SAXS_PATH].shape)

        # range-averaged ion-chamber monitors (upstream / downstream)
        _write_scalar(hf, INCIDENT_PATH, incident_avg)
        _write_scalar(hf, TRANSMITTED_PATH, transmitted_avg)

        if start_time is not None:
            hf[START_TIME][()] = start_time

        if FILE_LIST in hf:
            del hf[FILE_LIST]
        basenames = [os.path.basename(f) for f in included_files]
        hf.create_dataset(
            FILE_LIST,
            data=np.array(basenames, dtype=object),
            dtype=h5py.string_dtype(encoding='utf-8'),
        )


def process_group(group, file_ranges, timelist):
    """Average every frame range for a single group and write the output files."""
    flist_all = sorted(glob.glob(os.path.join(prefix, f'{group}*_results.hdf')))
    # ignore any previously produced averaged files (new 'Average_...' names as
    # well as older '..._Average_...' ones)
    flist_all = [f for f in flist_all if 'Average' not in os.path.basename(f)]
    assert flist_all, f'no dataset found in {prefix} for group {group}'
    frame_numbers = np.array([extract_frame(f) for f in flist_all])

    print(f'\n=== {group}: {len(flist_all)} raw files in {prefix} ===')
    for start, end in file_ranges:
        sel = (frame_numbers >= start) & (frame_numbers <= end)
        section_files = [f for f, m in zip(flist_all, sel) if m]
        assert section_files, f'no files found for {group} range ({start}, {end})'

        label = f'{group}_{start:05d}_{end:05d}'
        print(f'\nrange {start}-{end}: {len(section_files)} files')

        data_dict, read_files, incident, transmitted = read_range_data(section_files)
        assert read_files, f'no readable files for {group} range ({start}, {end})'

        mask = outlier_removal(data_dict, label=label, percentile=percentile)
        included_files = [f for f, keep in zip(read_files, mask) if keep]
        print(f'  kept {np.sum(mask)} of {len(mask)} files after outlier removal')

        avg_dict = average_datasets(data_dict=data_dict, mask=mask)

        # Average the two beam monitors over the SAME surviving files (the mask
        # from the g2/saxs outlier removal), so the averaged file's flux and
        # transmission are consistent with its averaged g2/saxs curves.
        incident_avg = float(np.nanmean(incident[mask]))
        transmitted_avg = float(np.nanmean(transmitted[mask]))

        # start time comes from the first (lowest frame number) file in the range
        start_time = get_start_time(section_files[0], timelist)
        if start_time is None:
            print(f'  WARNING: no time list entry for {os.path.basename(section_files[0])}')

        template = section_files[0]
        # replace the single frame number with the range, then prepend 'Average_'
        core = re.sub(
            r'_(\d+)_results\.hdf$',
            f'_{start:05d}_{end:05d}_results.hdf',
            os.path.basename(template),
        )
        out_name = f'Average_{core}'
        out_path = os.path.join(out_dir, out_name)

        save_average(template, out_path, avg_dict, start_time, included_files,
                     incident_avg, transmitted_avg)
        print(f'  saved -> {out_name}  (start_time={start_time}, '
              f'incident={incident_avg:.6g}, transmitted={transmitted_avg:.6g})')


def main():
    os.makedirs(out_dir, exist_ok=True)
    timelist = load_timelist(timelist_path)

    print(f'output dir: {out_dir}')
    for group, file_ranges in FILE_RANGES.items():
        process_group(group, file_ranges, timelist)


if __name__ == '__main__':
    main()

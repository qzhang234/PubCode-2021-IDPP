"""Reduce the raw 8-ID-I result files into the small set of files in ``data/``.

This is the ONE place where the beamline's reprocessed NeXus results are read.
Everything downstream -- saxpcs.py, saxs_evolution.py, g2_grid_SI.py,
thermal_cycle.py, contrast_calibration.py -- reads only what this script writes
into ``data/``, so the repository carries enough reduced data for the whole
analysis to be rerun without access to the beamline storage.

It produces three kinds of file:

1. AVERAGES over explicit frame-number ranges, one file per range, for every
   group listed in ``FILE_RANGES``.  This is the bulk of the output and the
   original purpose of the script.
2. The per-acquisition TEMPERATURE TRACE of the thermal-cycling sequence
   (``TRACE_GROUPS`` -> ``TRACE_CSV``).  Figure S6a plots the temperature of
   every one of its 2325 acquisitions, which is not an average and cannot be
   recovered from the averaged files.
3. A STACK of the individual repeats of the contrast standard
   (``STACK_GROUP`` -> ``STACK_NAME``).  Figure S4a shows all 50 repeats and
   their scatter, so an average would destroy the very thing it plots; the
   stack keeps each repeat's g2 on the raw NeXus paths, with the repeat index
   as the leading axis.

For each range this script:
  * globs the group's result files and selects those whose frame number is in
    the range,
  * removes outliers in two passes: a cross-correlation threshold on the SHAPE
    of g2, g2_err and saxs_1d, then a robust cut on acquisitions sitting far
    above the group median in the 0.004-0.008 A^-1 band (``spike_removal``),
    and finally any frames listed for that group in ``MANUAL_EXCLUDE``,
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
    as the originals -- a full copy of the first file in the range is used as
    the template, so every averaged file is a complete, self-describing NeXus
    file rather than a stripped extract.

The averaged file name starts with 'Average_' and carries the frame range so it
is easy to tell apart from the raw files, e.g.::

    Average_B0147_S3_7_300C10p_att00_Rq0_00951_01050_results.hdf

In the averaged file:
  * ``/entry/start_time`` is set to the acquisition time of the first file in
    the range (looked up in the beamline time list),
  * ``/xpcs/average/file_list`` is added, listing every file included in the
    average (i.e. the files that survived outlier removal).

Run it once, before any figure script::

    python average_ranges.py            # everything
    python average_ranges.py B0083      # just one group's averages
"""

import sys
import os
import re
import glob
import shutil

from datetime import datetime

import numpy as np
import h5py

sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..'))
from common.utils import outlier_removal, average_datasets

# --- PARAMETERS ---
# BEAMLINE-ONLY: the raw reprocessed NeXus results live here and nowhere else.
# This path is the reason average_ranges.py cannot run from a clone of the
# repository; every other 8-ID-I script reads only the data/ folder it writes.
prefix = '/home/8-id-i/2022-1/babnigg202203_nexus/reprocess_results'

# Groups to average: file-name header -> list of (start, end) frame ranges
# (inclusive), one averaged file per range.  Usually one range per group.

FILE_RANGES = {
    # --- Figures 3, S8, S9, S10: the isothermal series and its two references
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

    # --- Figure S6: seven thermal cycles of sample S1, 2022-03-04 -----------
    # Each cycle is a 6 C group of 50 acquisitions followed by a 270-acquisition
    # ramp to 34 C, of which two ten-acquisition windows near the top are used
    # (they land near 31.9 and 33.8 C).  A range simply selects whatever files
    # exist inside it, so the nominal (1, 50) is right even where a few
    # acquisitions are missing (B0078 has 49, D0077 has 39).
    'B0075': [(1, 50)], 'B0076': [(241, 250), (261, 270)],
    'B0078': [(1, 50)], 'B0079': [(241, 250), (261, 270)],
    'B0081': [(1, 50)], 'B0082': [(241, 250), (261, 270)],
    'B0083': [(1, 50)], 'B0084': [(241, 250), (261, 270)],
    'B0085': [(1, 50)], 'B0086': [(241, 250), (261, 270)],
    'B0087': [(1, 50)], 'B0088': [(241, 250), (261, 270)],
    'B0089': [(1, 50)], 'B0090': [(241, 250), (261, 270)],
    # the buffer for those cycles, held cold and measured twice
    'D0077': [(1, 50)],
    'D0080': [(1, 50)],
}

# --- SECOND-STAGE CUT: acquisitions far above the group at low q ---------
# outlier_removal() cuts on the SHAPE of log10 I(q) across the whole q range.
# An acquisition that is normal everywhere except the lowest q bins stays
# nearly parallel to the group mean and survives it.  That is exactly what a
# large object drifting through the 10 x 10 um beam during one 2 s exposure
# looks like: a low-q spike with ordinary counting statistics above 0.01 A^-1.
# Averaging such an acquisition into the group puts the scatterer, not the
# sample, into the group mean.  This second cut catches those.
#
# SPIKE_BAND is the q window the test is made in.  It is the same window
# Figure S6b already reports its intensities over, so the cut introduces no new
# choice of q range, and the per-acquisition mean of I(q) over that window is
# the test statistic.  The test is therefore on the intensity in one band, not
# on the low-q rise as such: an acquisition that is bright across the whole
# detector fails it too (B0147 frames 774 and 794, 2x the group median at low q
# and 4x at high q, are removed for that reason).  That is the intended
# behaviour -- neither kind of acquisition is measuring the sample -- but it is
# worth knowing that the band is where the test looks, not what it diagnoses.
#
# The test is a one-sided iterated modified z-score (Iglewicz & Hoaglin): an
# acquisition is dropped when its band mean lies more than SPIKE_Z robust
# standard deviations ABOVE the median, where the scale is 1.4826 x MAD and
# both the median and the MAD are recomputed from the survivors until the set
# stops changing.  Median and MAD are used rather than mean and standard
# deviation because a handful of 8x acquisitions inflate the ordinary
# statistics enough to hide themselves.  The test is one-sided because the
# artefact only ever ADDS scattering; a low reading means lost flux, which the
# monitor normalisation and the shape cut already handle.
#
# SPIKE_Z = 3 is the conventional value and is not tuned: in the worst group
# the same acquisitions are selected for any threshold from 2.5 to 3.0, and the
# selection is unchanged if SPIKE_BAND is moved to 0.0032-0.006 or 0.003-0.005.
#
# The cut removes 35 of the 1708 acquisitions that survive outlier_removal()
# (2.0 %).  It is nearly inert on the isothermal series behind Figures 3 and
# S8-S10 -- 2 acquisitions of B0147, which move that group's absolute-scale
# coefficient by 0.03 % and leave every fitted g2 parameter unchanged -- and
# does its work on the Figure S6 thermal cycles, where transient scatterers are
# common: 4, 2, 1, 11, 4, 1 and 1 acquisitions in the seven 6 C groups and 2
# and 7 in the two buffer groups, against none at all in any of the fourteen
# ten-acquisition 34 C windows.  See the "spike removal" note in
# thermal_cycle.py.
SPIKE_BAND = (0.004, 0.008)          # A^-1
SPIKE_Z = 3.0


def static_q(fname):
    """q of every column of saxs_1d, read from the file's own static q map.

    scattering_1d is stored as one value per (q, phi) bin, flattened; the index
    map gives each bin's position in the (dim0 = q, dim1 = phi) grid.
    """
    with h5py.File(fname, 'r') as hf:
        idx = hf['/xpcs/qmap/static_index_mapping'][()]
        q_list = hf['/xpcs/qmap/static_v_list_dim0'][()]
        n_phi = hf['/xpcs/qmap/static_v_list_dim1'].shape[0]
    return q_list[idx // n_phi]


def spike_removal(saxs_1d, q, mask, band=SPIKE_BAND, z=SPIKE_Z):
    """Flag acquisitions whose low-q band mean sits far above the group median.

    ``saxs_1d`` is (n_acquisitions, n_bins) as stored in the file, ``q`` the
    matching q of every bin from static_q(), and ``mask`` the acquisitions
    still standing after outlier_removal().  The statistic is built from the
    WHOLE group, not from ``mask``'s survivors, so the two cuts are independent
    tests of the same acquisitions rather than a chain in which the second one
    inherits the first one's percentile; ``mask`` is applied only at the end.
    Returns a keep-mask.
    """
    sel = (q >= band[0]) & (q < band[1])
    stat = np.nanmean(saxs_1d[:, sel], axis=1)
    keep = np.ones(len(stat), dtype=bool)
    for _ in range(20):
        med = np.median(stat[keep])
        scale = 1.4826 * np.median(np.abs(stat[keep] - med))
        if not scale > 0:
            break
        new = keep & (stat - med <= z * scale)
        if np.array_equal(new, keep):
            break
        keep = new
    return mask & keep


# Acquisitions dropped BY HAND, on top of the two automatic cuts:
#   group -> (start, end) -> [frame numbers to drop]
# Empty, and kept only as the escape hatch for a defect the two automatic cuts
# cannot see.  The eleven acquisitions of B0083 that were once listed here by
# eye are now removed by spike_removal(), which takes eleven acquisitions of
# that group, ten of them the same: it also takes frame 9, and frame 47 of the
# hand list survives.  Those two sit at the boundary of the group's own scatter
# and no ranking separates them, which is why the by-eye list was replaced.
MANUAL_EXCLUDE = {}

# Groups whose per-acquisition temperature and time are written to TRACE_CSV,
# in the order they were measured.  Figure S6a needs every acquisition, not the
# group averages, so this is an extraction rather than a reduction.  The
# temperature is that of the stage zone holding the SAMPLE (qnw1) for every
# group, buffer acquisitions included: panel (a) is the thermal history of the
# sample, and the buffer sat in a different zone that was held at 6 C
# throughout.
TRACE_GROUPS = ['B0075', 'B0076', 'D0077', 'B0078', 'B0079', 'D0080',
                'B0081', 'B0082', 'B0083', 'B0084', 'B0085', 'B0086',
                'B0087', 'B0088', 'B0089', 'B0090']
TRACE_CSV = 'thermal_cycle_temperature.csv'
TEMP_PATH = '/entry/sample/qnw1_temperature'

# The contrast standard: 50 repeat acquisitions of nano-porous glass.  Figure S4
# plots all 50, so they are stacked rather than averaged.
STACK_GROUP = 'F0145'
STACK_NAME = 'Stack_F0145_10nm_Glass_006C_att00_Rq0_00001_00050_results.hdf'
DELAY_PATH = '/xpcs/multitau/delay_list'
FRAME_TIME_PATH = '/entry/instrument/detector_1/frame_time'
DYN_Q_PATH = '/xpcs/qmap/dynamic_v_list_dim0'


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


def group_files(group):
    """Raw result files of a group, in frame order, ignoring earlier averages."""
    fl = sorted(glob.glob(os.path.join(prefix, f'{group}*_results.hdf')))
    fl = [f for f in fl if 'Average' not in os.path.basename(f)
          and 'Stack' not in os.path.basename(f)]
    assert fl, f'no dataset found in {prefix} for group {group}'
    return fl


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
    flist_all = group_files(group)
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
        print(f'  kept {np.sum(mask)} of {len(mask)} files after outlier removal')

        # transient low-q scatterers, invisible to the shape-based cut above
        spike = spike_removal(data_dict['saxs_1d'], static_q(read_files[0]), mask)
        n_spike = int(np.sum(mask & ~spike))
        if n_spike:
            hit = sorted(extract_frame(f) for f, m in zip(read_files, mask & ~spike) if m)
            print(f'  dropped {n_spike} more as low-q spikes {hit} '
                  f'-> {np.sum(spike)} remain')
        mask = spike

        # hand-listed frames, dropped on top of the automatic cuts
        drop = MANUAL_EXCLUDE.get(group, {}).get((start, end), [])
        if drop:
            hit = np.array([extract_frame(f) in drop for f in read_files])
            n_hit = int(np.sum(mask & hit))          # only those still standing
            mask = mask & ~hit
            print(f'  dropped {n_hit} more from the hand list {sorted(drop)} '
                  f'-> {np.sum(mask)} remain')

        included_files = [f for f, keep in zip(read_files, mask) if keep]

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


def write_temperature_trace(groups, timelist, out_path):
    """One row per acquisition: dataset, frame, elapsed seconds, temperature.

    The acquisition time comes from the time list rather than the file, because
    a 2025 reprocessing overwrote /entry/start_time with the reprocessing date.
    Elapsed time is measured from the first acquisition of the sequence.
    """
    rows = []
    for group in groups:
        for f in group_files(group):
            key = os.path.basename(f).replace('_results.hdf', '')
            stamp = timelist.get(key)
            if stamp is None:
                print(f'  WARNING: no time list entry for {key}')
                continue
            with h5py.File(f, 'r') as hf:
                T = float(np.asarray(hf[TEMP_PATH][()]).ravel()[0])
            rows.append([key, extract_frame(f),
                         datetime.strptime(stamp, '%Y-%m-%d %H:%M:%S'), T])
    rows.sort(key=lambda r: r[2])
    t0 = rows[0][2]
    with open(out_path, 'w') as fh:
        fh.write('dataset,frame,elapsed_s,qnw1_temperature_C\n')
        for key, frame, t, T in rows:
            fh.write(f'{key},{frame},{(t - t0).total_seconds():.0f},{T:.4f}\n')
    print(f'  {len(rows)} acquisitions, {(rows[-1][2] - t0).total_seconds() / 3600:.2f} h, '
          f'T = {min(r[3] for r in rows):.2f}-{max(r[3] for r in rows):.2f} C '
          f'-> {os.path.basename(out_path)}')


def write_stack(group, out_path):
    """Every repeat of a group in one file, repeat index as the leading axis.

    Only the fields Figure S4 uses are carried, on the same NeXus paths the raw
    files use, so the reading code is the same either way.
    """
    files = group_files(group)
    G, E = [], []
    for f in files:
        with h5py.File(f, 'r') as hf:
            G.append(hf[G2_PATH][()])
            E.append(hf[G2_ERR_PATH][()])
    with h5py.File(files[0], 'r') as hf, h5py.File(out_path, 'w') as out:
        for path in (DELAY_PATH, FRAME_TIME_PATH, DYN_Q_PATH):
            out.create_dataset(path, data=hf[path][()])
        out.create_dataset(G2_PATH, data=np.array(G))
        out.create_dataset(G2_ERR_PATH, data=np.array(E))
        out.create_dataset(FILE_LIST,
                           data=np.array([os.path.basename(f) for f in files],
                                         dtype=object),
                           dtype=h5py.string_dtype(encoding='utf-8'))
    print(f'  {len(files)} repeats, g2 {np.array(G).shape} '
          f'-> {os.path.basename(out_path)}')


def main(only=None):
    """Re-reduce every group, or only the ones named on the command line.

    ``python average_ranges.py`` rebuilds all of data/.  ``python
    average_ranges.py B0083`` rebuilds just that group's averages and leaves
    every other file, the temperature trace and the contrast stack untouched --
    useful when only one group's frame selection has changed.
    """
    os.makedirs(out_dir, exist_ok=True)
    timelist = load_timelist(timelist_path)

    print(f'output dir: {out_dir}')
    for group, file_ranges in FILE_RANGES.items():
        if only and group not in only:
            continue
        process_group(group, file_ranges, timelist)
    if only:
        print('\nsubset run: temperature trace and contrast stack left as they were')
        return

    print(f'\n=== temperature trace of the cycling sequence ===')
    write_temperature_trace(TRACE_GROUPS, timelist,
                            os.path.join(out_dir, TRACE_CSV))

    print(f'\n=== {STACK_GROUP}: individual repeats of the contrast standard ===')
    write_stack(STACK_GROUP, os.path.join(out_dir, STACK_NAME))


if __name__ == '__main__':
    main(only=set(sys.argv[1:]) or None)

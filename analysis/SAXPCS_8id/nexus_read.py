"""Readers for the averaged 8-ID-I NeXus files in data_2021-3/.

Every script that reads those files uses these functions, so every figure
reads them the same way.  The callers open the files read-only.
"""

import os
import re
from datetime import datetime

import numpy as np

# where each quantity lives inside an averaged file
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

# 'Average_H0404_D7_285C10p_att00_Rq0_00101_00200_results.hdf'
#     -> header 'H0404', first acquisition 101, last acquisition 200
NAME_PATTERN = re.compile(r'Average_([A-Za-z]\d+)_.*?_(\d+)_(\d+)_results')


def parse_name(fname):
    """(header, first acquisition, last acquisition) from an averaged file's name."""
    match = NAME_PATTERN.search(os.path.basename(fname))
    if match is None:
        raise ValueError(f'not the name of an averaged result file: {fname}')
    return match.group(1), int(match.group(2)), int(match.group(3))


def read_start_time(hf):
    """Start time of the first acquisition in the file, as a datetime."""
    text = hf[START_TIME_PATH][()].decode('utf-8')
    return datetime.strptime(text.strip(), TIME_FORMAT)


def read_iq(hf):
    """Time- and azimuthally averaged I(Q), in photons per pixel per frame.

    The file stores one intensity per (Q, phi) bin, flattened with phi running
    fastest, so integer division of the bin index by the number of phi sectors
    gives the Q index.  The entries sharing a Q index are averaged.
    """
    intensity = np.asarray(hf[SAXS_PATH][()], dtype=float).reshape(-1)
    bin_index = hf[STATIC_MAP_PATH][()]
    q_list = hf[STATIC_Q_PATH][()]
    n_phi = hf[STATIC_PHI_PATH].shape[0]
    q_index = bin_index // n_phi
    unique_q = np.unique(q_index)
    I = np.zeros(len(unique_q))
    for k in range(len(unique_q)):
        I[k] = np.nanmean(intensity[q_index == unique_q[k]])
    return q_list[unique_q], I


def read_g2(hf):
    """Correlation functions: (delay in s, g2, g2_err, Q of each dynamic bin).

    delay_list is stored in frames.  g2 and g2_err are indexed (delay, Q bin).
    """
    tau = hf[DELAY_PATH][()] * hf[FRAME_TIME_PATH][()]
    return tau, hf[G2_PATH][()], hf[G2_ERR_PATH][()], hf[DYN_Q_PATH][()]

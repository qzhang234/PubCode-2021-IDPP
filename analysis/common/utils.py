"""The first outlier cut and the group average of the 2021-3 reduction.

reduce_babnigg202110.py averages each group of acquisitions through
average_ranges.py, which is not part of this repository; it calls
outlier_removal() and average_datasets() from here.  Beamline-side: nothing on
the path from data_2021-3/ to a figure uses this module.
"""

import numpy as np


def keep_most_similar(curves, percentile=5, style='linear'):
    """Keep-mask that drops the acquisitions whose curve is least like the rest.

    curves has one row per acquisition.  Every pair of rows is compared by
    cosine similarity (dot product over the product of lengths: 1 for curves of
    identical shape, whatever their level).  Each acquisition is scored by the
    sum of its similarities, and the lowest `percentile` per cent are dropped.
    With style='log' the curves are compared as log10, zeros and negatives
    first set to 1.  NaNs are set to 0 on a copy, so the input is not changed.
    """
    x = np.array(curves, dtype=np.float64)
    x[np.isnan(x)] = 0
    if style == 'log':
        x[x <= 0] = 1
        x = np.log10(x)

    num = x.shape[0]
    similarity = np.zeros((num, num), dtype=np.float64)
    for n in range(num):
        for m in range(n + 1, num):
            value = np.sum(x[n] * x[m]) / np.sqrt(np.sum(x[n] ** 2) * np.sum(x[m] ** 2))
            similarity[n, m] = value
            similarity[m, n] = value

    score = np.sum(similarity, axis=1)
    return score >= np.percentile(score, percentile)


def outlier_removal(data_dict, label=None, percentile=5):
    """Keep-mask of the first outlier cut (SI section 5.3).

    data_dict holds 'g2', 'g2_err' and 'saxs_1d', each with one row per
    acquisition.  An acquisition is kept only if it passes keep_most_similar()
    on every one of them; saxs_1d is compared in log.  The cut sees the shape of
    a curve, not its level.  label names the group for the caller and is not
    used here.
    """
    keep = np.ones(len(data_dict['g2']), dtype=bool)
    for key in data_dict:
        if key == 'saxs_1d':
            style = 'log'
        else:
            style = 'linear'
        keep = keep & keep_most_similar(data_dict[key], percentile=percentile, style=style)
    return keep


def average_datasets(data_dict, mask):
    """Average the acquisitions that mask keeps, one curve per dataset.

    g2 and saxs_1d: the mean, ignoring NaNs.  g2_err: sqrt(sum sigma_i^2) / N,
    the uncertainty of the mean of N independent measurements, with N counted
    per point over the non-NaN entries.
    """
    average = {}
    for key in data_dict:
        kept = data_dict[key][mask]
        if key == 'g2_err':
            average[key] = np.sqrt(np.nansum(kept ** 2, axis=0))
            # divided in place, so the result keeps the float type of the input
            average[key] /= np.sum(~np.isnan(kept), axis=0)
        else:
            average[key] = np.nanmean(kept, axis=0)
    return average

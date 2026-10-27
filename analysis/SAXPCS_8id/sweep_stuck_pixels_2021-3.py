"""Find every stuck pixel in the 2021-3 (babnigg202110) detector data.

BEAMLINE-ONLY.  Needs /gdata, so it runs on adamite or amazonite, not kouga.

WHY IT EXISTS
-------------
The first pass, recorded in the docstring of make_qmap_babnigg202110.py, found
(30, 340) as "the only pixel on the detector above 50x its ring median".  It
missed (63, 166), which is stuck at about 1.5 photons per frame in every 2021-3
acquisition.  That pass used the two standards only and a 50x threshold.  This
one sweeps every scan prefix of the beamtime, including all four ELP loadings,
the empties, the air and the standards, and walks the threshold down until the
list of flagged pixels stops growing.

METHOD
------
xpcs/temporal_mean/scattering_2d in each result is the mean counts per frame at
each pixel, averaged over all 100,000 frames of that acquisition.  Averaging it
again over many acquisitions of a group beats the shot noise down far enough
that a pixel stuck at O(1) photons/frame stands out against a ring median of
O(1e-4).

A pixel is judged against the median of its own q ring, because the scattering
itself varies by orders of magnitude from the beamstop outward.  Rings are
fixed-width in q over currently unmasked pixels only.

A stuck pixel is one that is high in EVERY group, including the air and the
empty cells, where there is no sample to scatter.  A pixel high only in the ELP
data is sample scattering, not a detector fault.

Usage
-----
    python sweep_stuck_pixels_2021-3.py                 # sweep, report, write nothing
    python sweep_stuck_pixels_2021-3.py --nacq 12       # more acquisitions per group
"""

import os
import glob
import argparse
import collections

import numpy as np
import h5py

NEXUS = "/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus"
RESULTS = os.path.join(NEXUS, "reprocess_results")
QMAP = os.path.join(NEXUS, "babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm.hdf")

SHAPE = (512, 1024)
ENERGY = 10.945        # keV
DISTANCE = 7.8         # m
PIXEL_SIZE = 7.6e-5    # m

# Groups to sweep, as (label, scan prefix).  Every prefix of the beamtime is
# here: four ELP loadings at their hold temperatures and their 6 C references,
# the empties, the air, and the standards.
GROUPS = [
    ("ELP D8 28.00C", "H0402_D8_280C10p_att00_Rq0"),
    ("ELP D9 28.25C", "H0412_D9_282C10p_att00_Rq0"),
    ("ELP D6 28.25C", "H0407_D6_282C10p_att00_Rq0"),
    ("ELP D7 28.50C", "H0404_D7_285C10p_att00_Rq0"),
    ("ELP D9 28.75C", "H0413_D9_287C10p_att00_Rq0"),
    ("ELP D8 6C",     "H0403_D8_060C10p_att00_Rq0"),
    ("ELP D7 6C",     "H0404_D7_060C10p_att00_Rq0"),
    ("ELP D6 6C",     "H0407_D6_060C10p_att00_Rq0"),
    ("ELP D9 6C",     "H0413_D9_060C10p_att00_Rq0"),
    ("empty E0411",   "E0411_Empty_3_006C_att00_Rq0"),
    ("empty G0409",   "G0409_Empty_006C_att00_Rq0"),
    ("empty G0410",   "G0410_Empty_2_006C_att00_Rq0"),
    ("air E0408",     "E0408_Air_006C_att00_Rq0"),
    ("empty E0304",   "E0304_Empty_006C_att00_Rq0"),
    ("std 12IDB",     "B0002_12IDB_Standard_006C_att00_Rq0"),
    ("std silica",    "F2303_150Silica_Glycerol_006C_att00_Rq0"),
]

NO_SAMPLE = {"empty E0411", "empty G0409", "empty G0410", "air E0408", "empty E0304"}

# The 9-10 Oct calibration measurements.  They ran with the downstream ion
# chamber not recording and an I0 near 1e9 rather than 9e5, so the detector
# and the front end were set up differently from the 17 Oct ELP runs.
OCT_0910 = {"empty E0304", "std 12IDB", "std silica"}

# Already masked by the current qmap; carried through as a reference.
KNOWN_BAD = [(30, 340)]


def q_per_pixel(cx, cy):
    """|q| at every pixel, from the beamline geometry."""
    wavelength = 12.398419 / ENERGY                       # A
    y, x = np.indices(SHAPE)
    r = np.hypot((x - cx) * PIXEL_SIZE, (y - cy) * PIXEL_SIZE)
    two_theta = np.arctan2(r, DISTANCE)
    return 4.0 * np.pi / wavelength * np.sin(two_theta / 2.0)


def group_mean(prefix, nacq):
    """Mean counts/frame per pixel, averaged over the first nacq acquisitions."""
    files = sorted(glob.glob(os.path.join(RESULTS, prefix + "_*_results.hdf")))[:nacq]
    if not files:
        return None, 0
    total = np.zeros(SHAPE, dtype=np.float64)
    used = 0
    for path in files:
        try:
            with h5py.File(path, "r") as hf:
                total += np.asarray(hf["xpcs/temporal_mean/scattering_2d"][0], dtype=np.float64)
            used += 1
        except (OSError, KeyError):
            continue
    if not used:
        return None, 0
    return total / used, used


def ring_median(img, q, good, nring=400):
    """Median of each fixed-width q ring, broadcast back to every pixel."""
    qg = q[good]
    edges = np.linspace(qg.min(), qg.max(), nring + 1)
    idx = np.clip(np.digitize(q, edges) - 1, 0, nring - 1)
    med = np.full(nring, np.nan)
    for i in range(nring):
        sel = good & (idx == i)
        if sel.sum() >= 20:
            med[i] = np.median(img[sel])
    return med[idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--nacq", type=int, default=8,
                    help="acquisitions averaged per group")
    ap.add_argument("--thresholds", type=float, nargs="+",
                    default=[50, 20, 10, 5, 3],
                    help="ratio-to-ring-median thresholds to walk down")
    args = ap.parse_args()

    with h5py.File(QMAP, "r") as hf:
        mask = hf["qmap/mask"][()].astype(bool)
        cx = float(hf["qmap/beam_center_x"][()])
        cy = float(hf["qmap/beam_center_y"][()])
    print(f"qmap {os.path.basename(QMAP)}")
    print(f"  centre (x={cx}, y={cy}), {mask.sum()} good pixels "
          f"({mask.sum() / mask.size:.4f} valid fraction)\n")

    q = q_per_pixel(cx, cy)

    # --- per-group mean images and ratio maps --------------------------------
    images, ratios, nacq_used = {}, {}, {}
    print(f"{'group':16s} {'acq':>4s}  {'ring med':>10s}  {'max ratio':>10s}  at")
    for label, prefix in GROUPS:
        img, used = group_mean(prefix, args.nacq)
        if img is None:
            print(f"{label:16s} {'--':>4s}  (no results)")
            continue
        med = ring_median(img, q, mask)
        with np.errstate(invalid="ignore", divide="ignore"):
            ratio = np.where(mask & np.isfinite(med) & (med > 0), img / med, 0.0)
        images[label], ratios[label], nacq_used[label] = img, ratio, used
        iy, ix = np.unravel_index(np.argmax(ratio), SHAPE)
        print(f"{label:16s} {used:4d}  {med[iy, ix]:10.3e}  {ratio[iy, ix]:10.1f}  "
              f"({iy}, {ix})")

    if not ratios:
        raise SystemExit("no results read; is /gdata mounted?")

    ngroup = len(ratios)

    # --- walk the threshold down ---------------------------------------------
    # A ring median can be as small as one photon in the whole group, so a
    # ratio alone is not enough: also require the pixel to carry real counts.
    # MIN_CPF is well below the 1.5 photons/frame of the known pixel and well
    # above the 1e-4 of ordinary scattering out here.
    MIN_CPF = 0.01
    print(f"\nflagged pixels (ratio > threshold AND > {MIN_CPF} counts/frame), "
          f"by how many of the {ngroup} groups agree:")
    print(f"{'thresh':>7s}  {'>=all':>6s} {'>=3/4':>6s} {'>=half':>6s}  "
          f"pixels flagged in at least half")
    votes_at = {}
    for thr in sorted(args.thresholds, reverse=True):
        votes = np.zeros(SHAPE, dtype=int)
        for label in ratios:
            votes += (ratios[label] > thr) & (images[label] > MIN_CPF)
        votes_at[thr] = votes
        half = sorted(zip(*np.where(votes >= ngroup // 2)))
        shown = ", ".join(f"({y}, {x})" for y, x in half[:10])
        more = f" +{len(half) - 10} more" if len(half) > 10 else ""
        print(f"{thr:7.0f}  {int((votes >= ngroup).sum()):6d} "
              f"{int((votes >= 3 * ngroup // 4).sum()):6d} "
              f"{len(half):6d}  {shown}{more}")

    # Candidates: flagged in at least half the groups at the loosest threshold.
    thr_lo = min(args.thresholds)
    cands = sorted(zip(*np.where(votes_at[thr_lo] >= ngroup // 2)))
    # Always inspect the pixel the first pass already masks, as a reference.
    for ref in KNOWN_BAD:
        if ref not in cands:
            cands.append(ref)

    # --- per-group detail -----------------------------------------------------
    print("\nper-group detail for each candidate pixel "
          "(counts/frame, ring median, ratio):")
    for y, x in cands:
        tag = "  [already masked]" if not mask[y, x] else ""
        print(f"\n  pixel (y={y}, x={x})   q = {q[y, x]:.5f} A^-1{tag}")
        print(f"    {'group':16s} {'counts/frame':>13s} {'ring median':>13s} "
              f"{'ratio':>11s}")
        for label, _prefix in GROUPS:
            if label not in images:
                continue
            med = ring_median(images[label], q, mask)
            note = " (no sample)" if label in NO_SAMPLE else ""
            note += " (9-10 Oct)" if label in OCT_0910 else ""
            print(f"    {label:16s} {images[label][y, x]:13.4f} "
                  f"{med[y, x]:13.3e} {ratios[label][y, x]:11.1f}{note}")

    # --- neighbourhood, step 1 of the TODO -----------------------------------
    print("\n3x3 neighbourhood, mean counts/frame "
          "(centre pixel marked *):")
    for y, x in cands:
        print(f"\n  around (y={y}, x={x}):")
        for label in ("ELP D7 28.50C", "empty E0411", "air E0408", "std 12IDB"):
            if label not in images:
                continue
            img = images[label]
            print(f"    {label}")
            for dy in (-1, 0, 1):
                row = "      "
                for dx in (-1, 0, 1):
                    yy, xx = y + dy, x + dx
                    if 0 <= yy < SHAPE[0] and 0 <= xx < SHAPE[1]:
                        flag = "*" if (dy == 0 and dx == 0) else " "
                        row += f"{img[yy, xx]:11.5f}{flag}"
                    else:
                        row += f"{'--':>12s}"
                print(row)


if __name__ == "__main__":
    main()

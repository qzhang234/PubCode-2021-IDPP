"""Build the NeXus-era qmap for the 2021-3 isothermal gelation series.

BEAMLINE-ONLY.  Needs /gdata, which mounts on the 10.x hosts (amber) and not on
the 164.x analysis machines, so it cannot run from a clone of the repository.
It is kept here for provenance: it records exactly how the qmap was built.

WHICH QMAP TO USE
-----------------
Three files exist, in order of age:

    ..._Q0p0032_Lin_BadpixRm.hdf   what the default mode of this script writes
    ..._Q0p0032_Lin_CrossRm.hdf    built in the pySimpleMask 0.4.0 GUI; adds the
                                   beamstop cross arms and refines the beam
                                   centre to (x=510.5, y=259.5), 1.7 and 1.1
                                   pixels from the value recorded in the archive
    ..._Q0p0032_Lin_CrossRm2.hdf   CrossRm with the second stuck pixel masked

CrossRm2 is the one to process against.  `--crossrm2` derives it from CrossRm
by a surgical copy, for the reason given in make_crossrm2().  The default mode
is kept because it documents where the geometry, the mask sources and the
stuck pixels came from, all of which still hold.

WHY IT EXISTS
-------------
The 2021-3 beamtime (babnigg202110) holds the isothermal gelation series that
answers the referee question about the transition temperature: five samples
held at 28.00, 28.25, 28.50 and 28.75 C.  Those frames were never converted to
NeXus, so the current boost_corr, which reads NeXus only, cannot reach them.
The conversion needs a qmap, and the 2021-3 qmap that the Matlab pipeline used
(partitionMapLibrary/2021-3/babnigg202110_Rq0_S360_Sphi8_D18_Lin.h5) is gone.
This script rebuilds one in the modern format.

GEOMETRY
--------
Taken from the per-acquisition Data Exchange metadata of the 2021-3 archive:

    energy             10.945 keV
    detector distance  7.8 m
    pixel size         76 um
    detector           Rigaku XSPA-500k, 512 x 1024

The NeXus qmap indexes arrays as [y, x] with shape (ny=512, nx=1024), which is
the transpose of the Data Exchange convention, where beam_center_x runs along
the 512 axis.  The 2022-1 pair fixes the mapping: Data Exchange (x=235, y=584)
became qmap (x=581.5, y=235.5).  So the 2021-3 Data Exchange centre
(x=260.6, y=512.2) becomes qmap (x=512.2, y=260.6).

That centre is confirmed three ways.  It sits 1.6 px from the centroid of the
beamstop shadow in mask_beamstop_10_09.mat, so the direct beam is behind the
beamstop as it must be.  Feeding it through the geometry above reproduces the
smallest q the Matlab pipeline recorded, sqspan[0] = 0.00108568 A^-1, to six
figures.  And the beamstop shadow only lands on the beam when the mask is read
in this orientation, which fixes the transpose independently.

MASK
----
Three sources, unioned:

    beamstop        mask_beamstop_10_09.mat, the hand-drawn 2021-3 beamstop
                    (2548 px).  Read as-is: h5py already returns it [y, x].
    blemish         the current 8idRigaku500k bad-pixel map (52288 px): a
                    border frame, rows 0-18 and 494-511 and columns 0-14 and
                    1009-1023, plus 150 scattered bad pixels inside it.  The
                    detector has no module gaps and the map holds none: strip
                    the frame and not one interior row or column is fully
                    masked.  An earlier note here said it carried module gaps,
                    which was wrong and had been used to explain away a dip in
                    beta near 0.018 A^-1; that dip has no cause assigned.
    legacy bad px   the 47 pixels the 2021-3 Matlab mask excluded that the
                    current blemish does not, so nothing the original analysis
                    threw away comes back.
    stuck pixel     (30, 340), which none of the above catch.  This is the
                    BadpixRm part of the filename, the same step that removed
                    three pixels from the 2022-1 qmap.  Found by correlating
                    the two 2021-3 standards and looking for pixels far above
                    the median of their own q ring: it reads ~0.75 in the
                    12IDB standard, the silica and the ELP data alike, against
                    a ring median near 1e-4, while its immediate neighbours are
                    normal.  Left in, it alone drives the q = 0.0156 dynamic
                    bin to g2 ~ 2.4 at every delay in every sample.

    stuck pixel     (63, 166), q = 0.02143 A^-1, found 26 Oct 2026 and masked
                    in CrossRm2 only.  It reads 1.4982 to 1.5009 counts per
                    frame in every 17 Oct measurement alike: all four ELP
                    loadings at their hold temperatures, their 6 C references,
                    the three empty cells and the air.  Its ring median there
                    is 1.3e-6 in the ELP data, which is the one-photon
                    quantisation floor of an eight-acquisition average, and
                    8.8e-5 in the 12IDB standard.  Its eight neighbours read
                    1e-4 or less.  It lands in static bin Q = 0.0214 A^-1,
                    phi sector 7, where it spikes I(Q), and in dynamic ROI 19,
                    which it dominates.  The five fitted bins, 0.00359 to
                    0.00756 A^-1, are untouched, so f, D_med and the contrast
                    do not move.  It reads zero in 2022-1 (babnigg202203);
                    that beamtime is unaffected.

                    WHY THE FIRST PASS MISSED IT.  Not a loose threshold.  The
                    first pass correlated the two standards, B0002 and F2303,
                    and both were measured 9-10 Oct, when this pixel was still
                    healthy: it reads 0.0000 to 0.0001 counts per frame in
                    B0002, F2303 and E0304 alike, in line with its neighbours.
                    The fault appeared between 10 and 17 Oct, so no sweep over
                    the standards alone could have seen it.  Re-running the
                    sweep over every scan prefix of the beamtime, with the
                    threshold walked from 50x down to 3x its ring median,
                    returns (63, 166) and nothing else; see
                    sweep_stuck_pixels_2021-3.py.

then an inner cut at q = 0.0032 A^-1, matching the 2022-1 qmap so both datasets
start at the same q.  This geometry would allow 0.00109; the cut is a choice,
not a limit, and it costs ~8500 low-q pixels.

PARTITION
---------
Sq270 / Dq27 / Sphi8 / Dphi1, linear, matching the 2022-1 qmap rather than the
S360/Sphi8/D18 the Matlab pipeline used, so the two cycles bin the same way.

A NOTE ON REPRODUCIBILITY
-------------------------
The 2022-1 qmap was written by pysimplemask 0.2.0.post7, whose source is no
longer readable.  This script uses 0.3.1.post2.  Regenerating the 2022-1
partition from its own stored mask with 0.3.1 reproduces the index mappings
exactly and the q bin centres to 1.6e-5 relative: the two versions disagree
slightly in the q formula at the outer edge.  That is 0.0016% in q, far below
anything the measurement resolves, but it means this file is not bit-identical
to one 0.2.0.post7 would have written.

Usage
-----
    python make_qmap_babnigg202110.py                        # dry run, BadpixRm
    python make_qmap_babnigg202110.py --apply                # write BadpixRm
    python make_qmap_babnigg202110.py --crossrm2             # dry run, CrossRm2
    python make_qmap_babnigg202110.py --crossrm2 --apply     # write CrossRm2
"""

import os
import sys
import argparse

import types

import numpy as np
import h5py
from PIL import Image

PYSIMPLEMASK_SRC = "/home/beams10/8IDIUSER/Documents/Miaoqi/xpcs_tools/pySimpleMask_p2504/src"
sys.path.insert(0, PYSIMPLEMASK_SRC)

# pysimplemask/__init__.py imports the Qt GUI, which is not installed on the
# GPU hosts.  Only the two pure-numpy submodules below are needed, so register
# a stub package to import them without running that __init__.
if "pysimplemask" not in sys.modules:
    _pkg = types.ModuleType("pysimplemask")
    _pkg.__path__ = [os.path.join(PYSIMPLEMASK_SRC, "pysimplemask")]
    sys.modules["pysimplemask"] = _pkg

from pysimplemask.qmap import compute_transmission_qmap          # noqa: E402
from pysimplemask.utils import (                                  # noqa: E402
    generate_partition,
    combine_partitions,
    hash_numpy_dict,
    optimize_integer_array,
)

ARCHIVE = "/gdata/s8id-dmdtn/2021-3/babnigg202110"
BEAMSTOP_MAT = f"{ARCHIVE}/cluster_results/mask_beamstop_10_09.mat"
LEGACY_RESULT = (
    f"{ARCHIVE}/cluster_results/H0404_D7_285C10p_att00_Rq0/"
    "H0404_D7_285C10p_att00_Rq0_00001_0001-100000.hdf"
)
BLEMISH = "/home/beams/8IDIUSER/Documents/areaDetectorBlemish/8idRigaku500k/latest_blemish.tif"
SOURCE_BIN = (
    f"{ARCHIVE}/H0404_D7_285C10p_att00_Rq0/H0404_D7_285C10p_att00_Rq0_00001/"
    "H0404_D7_285C10p_att00_Rq0_00001.bin"
)

OUT_DIR = "/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus"
OUT_NAME = "babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_BadpixRm.hdf"

SHAPE = (512, 1024)          # (ny, nx), NeXus qmap orientation
BEAM_CENTER_X = 512.2        # along nx
BEAM_CENTER_Y = 260.6        # along ny
ENERGY = 10.945              # keV
DISTANCE = 7.8               # m
PIXEL_SIZE = 7.6e-5          # m
Q_INNER_CUT = 0.0032         # A^-1

# Stuck pixels this cycle, as (y, x) in NeXus orientation.  See the docstring.
EXTRA_BAD_PIXELS = [(30, 340), (63, 166)]

# The qmap actually processed against, and the one this script now derives.
CROSS_NAME = "babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm.hdf"
CROSS2_NAME = "babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm2.hdf"

# Pixels CrossRm2 masks that CrossRm does not.  Only the second stuck pixel:
# (30, 340) is already masked in CrossRm.
CROSS2_NEW_BAD = [(63, 166)]

SQ_NUM, SP_NUM = 270, 8
DQ_NUM, DP_NUM = 27, 1
STYLE = "linear"

# The q span the Matlab pipeline recorded for this cycle, used as a check.
# Agreement is ~2e-4 relative, which is a sub-0.03 pixel difference in the
# implied radius: the recorded centre and the one the Matlab qmap used differ
# by far less than a pixel.  The tolerance below is set to catch a wrong
# centre or a wrong transpose, both of which miss by orders of magnitude more.
SQSPAN_MATLAB = (0.0010856761, 0.0306956694)
SQSPAN_TOL = 5e-4


def build_mask(verbose=True):
    """Union the beamstop, the blemish and the legacy bad pixels."""
    with h5py.File(BEAMSTOP_MAT, "r") as f:
        usermask = f["usermask"][()]
    assert usermask.shape == SHAPE, f"beamstop mask is {usermask.shape}, expected {SHAPE}"
    beamstop_bad = usermask == 0

    blemish = np.asarray(Image.open(BLEMISH))
    assert blemish.shape == SHAPE
    blemish_bad = blemish == 0

    # xpcs/mask in the Matlab result is stored [x, y]; transpose to [y, x].
    with h5py.File(LEGACY_RESULT, "r") as f:
        legacy_bad = f["xpcs"]["mask"][()].T == 0
    assert legacy_bad.shape == SHAPE
    assert not (beamstop_bad & ~legacy_bad).any(), "beamstop not inside the legacy mask"
    legacy_only = legacy_bad & ~blemish_bad & ~beamstop_bad

    stuck = np.zeros(SHAPE, dtype=bool)
    for y, x in EXTRA_BAD_PIXELS:
        stuck[y, x] = True

    bad = beamstop_bad | blemish_bad | legacy_bad | stuck
    if verbose:
        print(f"  beamstop                 {beamstop_bad.sum():7d}")
        print(f"  blemish                  {blemish_bad.sum():7d}")
        print(f"  legacy-only bad pixels   {legacy_only.sum():7d}")
        print(f"  stuck pixels             {int((stuck & ~blemish_bad).sum()):7d}")
        print(f"  union before q cut       {bad.sum():7d}")
    return ~bad



def make_crossrm2(out_dir, apply_changes, overwrite):
    """Derive CrossRm2 from CrossRm by masking the pixels in CROSS2_NEW_BAD.

    WHY A SURGICAL COPY RATHER THAN A REBUILD
    -----------------------------------------
    CrossRm was written by pysimplemask 0.4.0.post2.  Rebuilding its partitions
    from its own stored mask with the version installed here reproduces the
    index mappings exactly and the q bin centres to 1.6e-5 relative, the same
    q-formula difference between versions that the docstring records for the
    2022-1 qmap.  That tiny shift moves 566 of 453,701 static pixels and 34
    dynamic pixels across a bin edge.  Nothing the measurement resolves, but it
    is a change beyond the one intended.  Copying CrossRm and clearing the
    stuck pixel changes exactly the pixel in question and leaves every q value,
    every bin edge and every index mapping byte-identical.

    THE INVARIANT
    -------------
    In CrossRm, a masked pixel has mask False and both roi maps 0, and the
    counts agree exactly: 70,587 masked pixels, 70,587 zeros in each roi map.
    The blemish dataset is the raw detector TIFF and is NOT touched; (30, 340)
    has blemish 1 there.  CrossRm2 follows the same convention.
    """
    src = os.path.join(out_dir, CROSS_NAME)
    dst = os.path.join(out_dir, CROSS2_NAME)
    if not os.path.exists(src):
        raise SystemExit(f"missing source qmap {src}")

    def decode(v):
        """h5py hands back bytes for string datasets; hash_numpy_dict needs str."""
        if isinstance(v, bytes):
            return v.decode()
        if isinstance(v, np.ndarray) and v.dtype == object:
            return [e.decode() if isinstance(e, bytes) else e for e in v]
        return v

    with h5py.File(src, "r") as hf:
        grp = hf["qmap"]
        data = {k: decode(grp[k][()]) for k in grp.keys()}
        attrs = dict(grp.attrs)
        units = {k: dict(grp[k].attrs) for k in grp.keys() if "_v_list_dim" in k}

    mask = data["mask"].astype(bool)
    srm, drm = data["static_roi_map"], data["dynamic_roi_map"]
    n_masked0 = int((~mask).sum())
    assert int((srm == 0).sum()) == n_masked0 == int((drm == 0).sum()), \
        "CrossRm mask and roi maps disagree; refusing to edit"

    print(f"source  {CROSS_NAME}")
    print(f"  good pixels {int(mask.sum())}, masked {n_masked0}, "
          f"valid fraction {mask.sum() / mask.size:.4f}")

    changed = []
    for y, x in CROSS2_NEW_BAD:
        if not mask[y, x]:
            print(f"  ({y}, {x}) already masked; nothing to do")
            continue
        sroi, droi = int(srm[y, x]), int(drm[y, x])
        n_s = int((srm == sroi).sum())
        n_d = int((drm == droi).sum())
        print(f"  masking ({y}, {x}): static ROI {sroi} loses 1 of {n_s} px, "
              f"dynamic ROI {droi} loses 1 of {n_d} px")
        assert n_s > 1 and n_d > 1, "masking would empty a ROI"
        mask[y, x] = False
        srm[y, x] = 0
        drm[y, x] = 0
        changed.append((y, x))

    data["mask"] = mask
    data["static_roi_map"] = srm
    data["dynamic_roi_map"] = drm

    n_masked1 = int((~mask).sum())
    assert int((srm == 0).sum()) == n_masked1 == int((drm == 0).sum())
    assert n_masked1 == n_masked0 + len(changed)
    frac = mask.sum() / mask.size
    print(f"result  {CROSS2_NAME}")
    print(f"  good pixels {int(mask.sum())}, masked {n_masked1}, "
          f"valid fraction {frac:.4f}")

    # Step 4 of the reprocessing TODO: the runner passes
    # --crop-ratio-threshold 0.5, so the valid fraction must stay above 0.5.
    assert frac > 0.5, "valid fraction fell below the crop threshold"
    print(f"  valid fraction {frac:.4f} > 0.50 crop threshold: OK")

    # Every ROI that CrossRm kept must still be populated, or the index
    # mappings, which are copied verbatim, would point at an empty bin.  The
    # roi maps hold compacted ids 1..N with 0 for masked, while the index
    # mapping holds the id each one had in the full 270x8 (or 27x1) grid, so
    # the check is against 1..N and not against the mapping values.
    for name, roi_map, mapping in (("static", srm, data["static_index_mapping"]),
                                   ("dynamic", drm, data["dynamic_index_mapping"])):
        n_roi = np.asarray(mapping).size
        present = np.unique(roi_map)
        present = present[present != 0]
        missing = sorted(set(range(1, n_roi + 1)) - set(int(i) for i in present))
        assert not missing, f"{name} ROIs now empty: {missing}"
        print(f"  all {n_roi} {name} ROIs still populated")

    if not apply_changes:
        print(f"\ndry run; would write {dst}")
        return

    for key, val in data.items():
        data[key] = optimize_integer_array(val)
    hash_val = hash_numpy_dict(data)

    if os.path.exists(dst) and not overwrite:
        raise SystemExit(f"refusing to overwrite existing {dst}; pass --overwrite")

    tmp = dst + ".tmp"
    with h5py.File(tmp, "w") as hf:
        grp = hf.create_group("/qmap")
        for key, val in data.items():
            compression = "lzf" if isinstance(val, np.ndarray) and val.size > 1024 else None
            dset = grp.create_dataset(key, data=val, compression=compression)
            for ak, av in units.get(key, {}).items():
                dset.attrs[ak] = av
        for ak, av in attrs.items():
            grp.attrs[ak] = av
        grp.attrs["hash"] = hash_val
    os.replace(tmp, dst)

    print(f"\nwrote {dst}")
    print(f"hash {hash_val}  (CrossRm was {attrs.get('hash')})")
    print(f"pysimplemask version attr carried over: {attrs.get('version')}")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write the qmap")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--overwrite", action="store_true",
                    help="replace an existing qmap of the same name")
    ap.add_argument("--crossrm2", action="store_true",
                    help="derive CrossRm2 from CrossRm by masking CROSS2_NEW_BAD; "
                         "this is the qmap to process against")
    args = ap.parse_args()

    if args.crossrm2:
        make_crossrm2(args.out_dir, args.apply, args.overwrite)
        return

    print("mask components (pixels excluded):")
    mask = build_mask()

    qmap, qmap_unit = compute_transmission_qmap(
        ENERGY, (BEAM_CENTER_Y, BEAM_CENTER_X), SHAPE, PIXEL_SIZE, DISTANCE, 0.0
    )
    q = qmap["q"]

    # Check the geometry against the Matlab pipeline over the mask the Matlab
    # pipeline actually used, which is the legacy mask, not the one built above.
    with h5py.File(LEGACY_RESULT, "r") as f:
        legacy_good = f["xpcs"]["mask"][()].T != 0
    print("\ngeometry check against the Matlab qmap, over the legacy mask:")
    for label, got, want in (
        ("min q", q[legacy_good].min(), SQSPAN_MATLAB[0]),
        ("max q", q[legacy_good].max(), SQSPAN_MATLAB[1]),
    ):
        rel = abs(got - want) / want
        print(f"  {label}: rebuilt {got:.8f}  Matlab {want:.8f}  rel {rel:.2e}")
        assert rel < SQSPAN_TOL, "beam centre or geometry disagrees with the Matlab pipeline"

    mask = mask & (q >= Q_INNER_CUT)
    print(f"\nafter the q >= {Q_INNER_CUT} cut: {mask.sum()} good pixels, {(~mask).sum()} excluded")
    print(f"q range retained: {q[mask].min():.6f} .. {q[mask].max():.6f} A^-1")

    sq = generate_partition("q", mask, q, SQ_NUM, style=STYLE, phi_offset=None)
    sp = generate_partition("phi", mask, qmap["phi"], SP_NUM, style=STYLE,
                            phi_offset=0.0, symmetry_fold=1)
    dq = generate_partition("q", mask, q, DQ_NUM, style=STYLE, phi_offset=None)
    dp = generate_partition("phi", mask, qmap["phi"], DP_NUM, style=STYLE,
                            phi_offset=0.0, symmetry_fold=1)
    static_map = combine_partitions(sq, sp, prefix="static")
    dynamic_map = combine_partitions(dq, dp, prefix="dynamic")

    partition = {
        "beam_center_x": BEAM_CENTER_X,
        "beam_center_y": BEAM_CENTER_Y,
        "pixel_size": PIXEL_SIZE,
        "mask": mask,
        "energy": ENERGY,
        "detector_distance": DISTANCE,
        "map_names": ["q", "phi"],
        "map_units": [qmap_unit["q"], qmap_unit["phi"]],
        "source_file": SOURCE_BIN,
    }
    partition.update(static_map)
    partition.update(dynamic_map)

    print(f"\nstatic  q bins {SQ_NUM} x phi {SP_NUM}: "
          f"{sq['v_list'][0]:.6f} .. {sq['v_list'][-1]:.6f} A^-1")
    print(f"dynamic q bins {DQ_NUM} x phi {DP_NUM}: "
          f"{dq['v_list'][0]:.6f} .. {dq['v_list'][-1]:.6f} A^-1")
    print(f"static ROIs kept  {static_map['static_index_mapping'].size}")
    print(f"dynamic ROIs kept {dynamic_map['dynamic_index_mapping'].size}")

    out_path = os.path.join(args.out_dir, OUT_NAME)
    if not args.apply:
        print(f"\ndry run; would write {out_path}")
        return

    for key, val in partition.items():
        partition[key] = optimize_integer_array(val)
    hash_val = hash_numpy_dict(partition)

    os.makedirs(args.out_dir, exist_ok=True)
    if os.path.exists(out_path) and not args.overwrite:
        raise SystemExit(f"refusing to overwrite existing {out_path}; pass --overwrite")

    with h5py.File(out_path, "w") as hf:
        grp = hf.create_group("/qmap")
        for key, val in partition.items():
            compression = "lzf" if isinstance(val, np.ndarray) and val.size > 1024 else None
            dset = grp.create_dataset(key, data=val, compression=compression)
            if "_v_list_dim" in key:
                dim = int(key[-1])
                dset.attrs["unit"] = partition["map_units"][dim]
                dset.attrs["name"] = partition["map_names"][dim]
                dset.attrs["size"] = np.asarray(val).size
        grp.attrs["hash"] = hash_val
        grp.attrs["version"] = "0.3.1.post2"

    print(f"\nwrote {out_path}")
    print(f"hash {hash_val}")


if __name__ == "__main__":
    main()

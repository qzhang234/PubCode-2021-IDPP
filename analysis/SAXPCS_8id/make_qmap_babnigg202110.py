"""Build the NeXus-era qmap (q partition and mask) of the 2021-3 isothermal series.

BEAMLINE-ONLY: needs /gdata (amber, adamite, amazonite) and the pySimpleMask
source.  Kept for provenance: it records how the qmap was built.

    python make_qmap_babnigg202110.py                        # dry run, BadpixRm
    python make_qmap_babnigg202110.py --apply                # write BadpixRm
    python make_qmap_babnigg202110.py --crossrm2             # dry run, CrossRm2
    python make_qmap_babnigg202110.py --crossrm2 --apply     # write CrossRm2

The series is correlated against ..._Lin_CrossRm2.hdf.  The default mode
writes ..._Lin_BadpixRm.hdf from the geometry and the mask sources below;
CrossRm (built in the pySimpleMask 0.4.0 GUI, which adds the beamstop cross
arms and refines the beam centre) is not reproduced by this script, and
--crossrm2 derives CrossRm2 from it by masking the second stuck pixel only
(make_crossrm2).  Geometry: 10.945 keV, 7.8 m, 76 um pixels, 512 x 1024, beam
centre (x=512.2, y=260.6) in NeXus orientation, the transpose of the Data
Exchange value.  Mask: the 2021-3 beamstop, the detector blemish map, the
pixels the 2021 Matlab mask excluded, and the stuck pixels (30, 340) and
(63, 166); then q >= 0.0032 A^-1.  Partition: 270 x 8 static, 27 x 1 dynamic,
linear.  The repository README gives the provenance in more detail.
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

# pysimplemask/__init__.py imports the Qt GUI, which the GPU hosts do not
# have.  Only two pure-numpy submodules are needed, so an empty package object
# pointing at the source directory is registered first; importing the
# submodules then skips that __init__.  (Unavoidable without the GUI.)
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


def compression_for(val):
    """lzf for arrays of more than 1024 entries, none otherwise."""
    if isinstance(val, np.ndarray) and val.size > 1024:
        return "lzf"
    return None


def build_mask():
    """Union of the beamstop, the blemish, the legacy bad pixels and the stuck
    pixels; True where a pixel is kept.  Prints the size of each part."""
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
    print(f"  beamstop                 {beamstop_bad.sum():7d}")
    print(f"  blemish                  {blemish_bad.sum():7d}")
    print(f"  legacy-only bad pixels   {legacy_only.sum():7d}")
    print(f"  stuck pixels             {int((stuck & ~blemish_bad).sum()):7d}")
    print(f"  union before q cut       {bad.sum():7d}")
    return ~bad

def decode_strings(v):
    """h5py returns bytes for string datasets; hash_numpy_dict needs str."""
    if isinstance(v, bytes):
        return v.decode()
    if isinstance(v, np.ndarray) and v.dtype == object:
        decoded = []
        for e in v:
            if isinstance(e, bytes):
                decoded.append(e.decode())
            else:
                decoded.append(e)
        return decoded
    return v


def make_crossrm2(out_dir, apply_changes, overwrite):
    """Derive CrossRm2 from CrossRm by masking the pixels in CROSS2_NEW_BAD.

    A copy rather than a rebuild: rebuilding CrossRm's partitions with the
    installed pySimpleMask shifts the q bin centres by 1.6e-5 relative (a
    difference between versions) and moves 566 static and 34 dynamic pixels
    across a bin edge.  The copy changes only the pixel in question.  In
    CrossRm a masked pixel has mask False and both roi maps 0 (70,587 of
    each); the blemish dataset is the raw detector map and is not touched.
    """
    src = os.path.join(out_dir, CROSS_NAME)
    dst = os.path.join(out_dir, CROSS2_NAME)
    if not os.path.exists(src):
        raise SystemExit(f"missing source qmap {src}")

    data = {}
    units = {}
    with h5py.File(src, "r") as hf:
        grp = hf["qmap"]
        for k in grp.keys():
            data[k] = decode_strings(grp[k][()])
            if "_v_list_dim" in k:
                units[k] = dict(grp[k].attrs)
        attrs = dict(grp.attrs)

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
        present = set()
        for i in np.unique(roi_map):
            present.add(int(i))
        missing = []
        for roi in range(1, n_roi + 1):
            if roi not in present:
                missing.append(roi)
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
            dset = grp.create_dataset(key, data=val, compression=compression_for(val))
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
            dset = grp.create_dataset(key, data=val, compression=compression_for(val))
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

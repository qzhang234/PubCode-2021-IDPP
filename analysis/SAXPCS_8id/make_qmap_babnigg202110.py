"""Build the NeXus-era qmap for the 2021-3 isothermal gelation series.

BEAMLINE-ONLY.  Needs /gdata, which mounts on the 10.x hosts (amber) and not on
the 164.x analysis machines, so it cannot run from a clone of the repository.
It is kept here for provenance: it records exactly how the qmap was built.

SUPERSEDED
----------
The qmap this script writes, ..._Q0p0032_Lin_BadpixRm.hdf, was replaced by
..._Q0p0032_Lin_CrossRm.hdf, built in the pySimpleMask 0.4.0 GUI.  That one
also masks the beamstop cross arms and refines the beam centre to
(x=510.5, y=259.5), 1.7 and 1.1 pixels from the value recorded in the
archive.  The CrossRm file is the one to process against.  This script is
kept because it documents where the geometry, the mask sources and the
stuck pixel came from, all of which still hold.

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
    blemish         the current 8idRigaku500k bad-pixel map (52288 px),
                    which carries the module gaps and the detector border.
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
                    normal.  It is the only pixel on the detector above 50x its
                    ring median.  Left in, it alone drives the q = 0.0156
                    dynamic bin to g2 ~ 2.4 at every delay in every sample.

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
    python make_qmap_babnigg202110.py            # dry run; reports, writes nothing
    python make_qmap_babnigg202110.py --apply    # write the qmap
"""

import os
import sys
import argparse

import numpy as np
import h5py
import tifffile

PYSIMPLEMASK_SRC = "/home/beams10/8IDIUSER/Documents/Miaoqi/xpcs_tools/pySimpleMask_p2504/src"
sys.path.insert(0, PYSIMPLEMASK_SRC)
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
EXTRA_BAD_PIXELS = [(30, 340)]

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

    blemish = np.asarray(tifffile.imread(BLEMISH))
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


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true", help="write the qmap")
    ap.add_argument("--out-dir", default=OUT_DIR)
    ap.add_argument("--overwrite", action="store_true",
                    help="replace an existing qmap of the same name")
    args = ap.parse_args()

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

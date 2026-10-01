"""Convert the 2021-3 isothermal gelation series from Data Exchange to NeXus.

BEAMLINE-ONLY.  Needs /gdata, which mounts on the 10.x hosts (amber) and not on
the 164.x analysis machines, so it cannot run from a clone of the repository.
It is kept here for provenance: it records exactly how the NeXus tree was made.

WHAT IT DOES
------------
For every acquisition under the H04* scan groups of babnigg202110 it writes

    babnigg202110_nexus/<acq>/<acq>.bin             symlink to the archive frame
    babnigg202110_nexus/<acq>/<acq>_metadata.hdf    NeXus metadata

mirroring the layout of babnigg202203_nexus, which is what the current
boost_corr expects.  The archive is opened read-only and never modified.

HOW IT DIFFERS FROM THE 2022-1 CONVERSION
-----------------------------------------
The 2025 rewrite that produced babnigg202203_nexus left template constants in
four places that were actually measured, and restore_dx_metadata.py had to go
back afterwards and copy the originals across.  This script reads those values
from the archive in the first place, using the same mapping that script
established:

    incident_beam_intensity     <- source_end/I0Monitor
    transmitted_beam_intensity  <- source_end/TransmissionMonitor
    ring_current                <- source_end/current
    start_time, end_time        <- source_begin/datetime, source_end/datetime

so no repair pass is needed.  It also writes the real geometry, detector name
and sample temperatures rather than leaving the template values in place.

A NOTE ON THE BEAM CENTRE
-------------------------
Data Exchange indexes the detector with beam_center_x along the 512 axis;
NeXus indexes it [y, x] with x along the 1024 axis.  The two are transposes of
each other, so the archive's (x=260.6, y=512.2) becomes (x=512.2, y=260.6).
That value was then refined in pySimpleMask against the measured scattering to
(x=510.5, y=259.5), which is what this script writes and what the delivered
qmap uses.

The pair beam_center_position_x/y is NOT a second copy of the beam centre.
The APS_8IDI reader treats it as ccdx0/ccdy0, the detector position at which
the beam centre was calibrated, and reconstructs the centre as

    beam_center_x = beam_center_x + (position_x - beam_center_position_x)
                                    / x_pixel_size

So this script writes beam_center_position_x = position_x and
beam_center_position_y = position_y.  The correction is then exactly zero and
the centre is read back as written.  Setting them equal also sidesteps a unit
mismatch in that expression: the Data Exchange stage positions are millimetres
while x_pixel_size is metres, so any non-zero difference would be scaled by
1000.  The archive's own stage_x and stage_zero_x differ by 0.0004 mm, a
twentieth of a pixel, so nothing measured is lost by setting them equal.

Usage
-----
    python nexus_convert_babnigg202110.py                 # dry run
    python nexus_convert_babnigg202110.py --apply
    python nexus_convert_babnigg202110.py --apply --limit 20     # smoke test
"""

import os
import re
import csv
import sys
import glob
import copy
import uuid
import argparse
import multiprocessing
from datetime import datetime

import numpy as np
import h5py

from nexus_xpcs_aps.nexus_schema import xpcs_schema
from nexus_xpcs_aps.nexus_utils import update_schema_at_runtime, create_nexus_format_metadata

ARCHIVE = "/gdata/s8id-dmdtn/2021-3/babnigg202110"
DEST = "/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus"
SCAN_GLOB = "H04*"
MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "nexus_convert_babnigg202110_manifest.csv")

DX_TIME_FORMAT = "%a %b %d %H:%M:%S %Y"      # 'Sun Oct 17 15:38:32 2021'
OUT_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"        # what nexus_read.py parses

# Geometry, in NeXus orientation and SI units.  See the module docstring.
# Refined in pySimpleMask against the measured scattering.
BEAM_CENTER_X, BEAM_CENTER_Y = 510.5, 259.5
DETECTOR_NAME = "rigaku500k"
CYCLE, PROPOSAL = "2021-3", "babnigg202110"

# Sample label and the QNW zone that holds it, read off the scan name.
SCAN_RE = re.compile(r"^(H\d{4})_(D\d+)_(\d{3})C10p_att00_Rq0$")


def scalar(obj):
    return float(np.asarray(obj[()]).ravel()[0])


def optional(group, name, missing):
    """A measured scalar, or NaN if this acquisition did not record it.

    The Data Exchange schema drifted during the 2021-3 beamtime: the Oct 9
    scans carry no source_end/TransmissionMonitor, the Oct 17 ones do.  NaN
    rather than the schema's placeholder constant, so a quantity that was
    never measured cannot be mistaken downstream for one that was.
    """
    if name in group:
        return scalar(group[name])
    missing.append(name)
    return float("nan")


def text(obj):
    v = np.asarray(obj[()]).ravel()[0]
    return v.decode("utf-8") if isinstance(v, bytes) else str(v)


def dx_file(acq_dir, acq):
    """The Data Exchange metadata file of one acquisition, or None."""
    direct = os.path.join(acq_dir, f"{acq}_0001-100000.hdf")
    if os.path.exists(direct):
        return direct
    found = sorted(glob.glob(os.path.join(acq_dir, "*_0001-*.hdf")))
    return found[0] if found else None


def read_metadata(path, acq):
    """Map one acquisition's Data Exchange record onto NeXus paths.

    Returns (metadata, missing), where missing names the measured fields this
    acquisition did not record.  A missing field is NaN, never a constant.
    """
    missing = []
    with h5py.File(path, "r") as hf:
        acqg = hf["measurement/instrument/acquisition"]
        det = hf["measurement/instrument/detector"]
        beg = hf["measurement/instrument/source_begin"]
        end = hf["measurement/instrument/source_end"]
        smp = hf["measurement/sample"]

        def stamp(obj):
            return datetime.strptime(text(obj).strip(), DX_TIME_FORMAT).strftime(OUT_TIME_FORMAT)

        translation = np.asarray(smp["translation"][()]).ravel()
        md = {
            "/entry/entry_identifier": acq,
            "/entry/entry_identifier_uuid": uuid.uuid4().hex,
            "/entry/beamline": "APS-8-ID-I",
            "/entry/scan_number": int(scalar(acqg["specscan_data_number"])),
            "/entry/start_time": stamp(beg["datetime"]),
            "/entry/end_time": stamp(end["datetime"]),
            "/entry/user/cycle": CYCLE,
            "/entry/user/proposal_id": PROPOSAL,
            # The Data Exchange archive does not record these.  "unknown" is
            # written rather than the schema's "John Doe" placeholder, which
            # would read as real.
            "/entry/user/name": "unknown",
            "/entry/user/email": "unknown",
            "/entry/user/institution": "unknown",

            "/entry/instrument/detector_1/detector_name": DETECTOR_NAME,
            "/entry/instrument/detector_1/count_time": scalar(det["exposure_time"]),
            "/entry/instrument/detector_1/frame_time": scalar(det["exposure_period"]),
            # Data Exchange records these in mm; NeXus wants m.
            "/entry/instrument/detector_1/distance": scalar(det["distance"]) * 1e-3,
            "/entry/instrument/detector_1/x_pixel_size": scalar(det["x_pixel_size"]) * 1e-3,
            "/entry/instrument/detector_1/y_pixel_size": scalar(det["y_pixel_size"]) * 1e-3,
            "/entry/instrument/detector_1/beam_center_x": BEAM_CENTER_X,
            "/entry/instrument/detector_1/beam_center_y": BEAM_CENTER_Y,
            # ccdx0/ccdy0 equal ccdx/ccdy, so the reader adds no offset.
            "/entry/instrument/detector_1/position_x": scalar(acqg["stage_x"]),
            "/entry/instrument/detector_1/position_y": scalar(acqg["stage_z"]),
            "/entry/instrument/detector_1/beam_center_position_x": scalar(acqg["stage_x"]),
            "/entry/instrument/detector_1/beam_center_position_y": scalar(acqg["stage_z"]),
            "/entry/instrument/detector_1/compression": "none",

            "/entry/instrument/incident_beam/incident_energy": scalar(beg["energy"]),
            "/entry/instrument/incident_beam/incident_beam_intensity": optional(end, "I0Monitor", missing),
            "/entry/instrument/incident_beam/transmitted_beam_intensity": optional(end, "TransmissionMonitor", missing),
            "/entry/instrument/incident_beam/ring_current": optional(end, "current", missing),

            "/entry/instrument/attenuator_1/attenuator_transmission": scalar(acqg["attenuation"]),
            "/entry/instrument/bluesky/spec_file": text(acqg["specfile"]),
            "/entry/instrument/bluesky/parent_folder": text(acqg["root_folder"]),
            "/entry/instrument/bluesky/bluesky_plan": "none",
            "/entry/instrument/bluesky/bluesky_plan_kwargs": "none",
            "/entry/instrument/bluesky/bluesky_version": "none",
            "/entry/instrument/bluesky/scan_id": int(scalar(acqg["specscan_data_number"])),
            "/entry/instrument/datamanagement/workflow_name": "nexus_convert_babnigg202110.py",
            "/entry/instrument/datamanagement/workflow_version": "1",
            "/entry/instrument/datamanagement/workflow_kwargs": "{}",

            "/entry/sample/position_x": float(translation[0]),
            "/entry/sample/position_y": float(translation[1]),
            "/entry/sample/position_z": float(translation[2]),
            "/entry/sample/qnw_lakeshore": optional(smp, "temperature_A", missing),
        }
        for z in (1, 2, 3):
            md[f"/entry/sample/qnw{z}_temperature"] = optional(smp, f"QNW_Zone{z}_Temperature", missing)
            md[f"/entry/sample/qnw{z}_temperature_set"] = optional(smp, f"QNW_Zone{z}_Temperature_Set", missing)
    return md, missing


def describe(scan):
    """Short and full sample descriptions from the scan name."""
    m = SCAN_RE.match(scan)
    if not m:
        return scan, scan
    _, sample, temp = m.groups()
    # The three digits are the Zone 3 setpoint times ten, truncated:
    # 282 is 28.25 C, 287 is 28.75 C.  The readback in the file is authoritative.
    return sample, f"ELP sample {sample}, nominal hold {int(temp) / 10.0:.2f} C"


def convert(job):
    """Convert one acquisition.  Returns a manifest row, or an error row."""
    scan, acq, acq_dir = job
    src_bin = os.path.join(acq_dir, f"{acq}.bin")
    out_dir = os.path.join(DEST, acq)
    out_bin = os.path.join(out_dir, f"{acq}.bin")
    out_md = os.path.join(out_dir, f"{acq}_metadata.hdf")

    dx = dx_file(acq_dir, acq)
    if dx is None:
        return (scan, acq, "", "", "", "", "no Data Exchange metadata")
    if not os.path.exists(src_bin):
        return (scan, acq, "", "", "", "", "no .bin in the archive")

    try:
        md, missing = read_metadata(dx, acq)
    except Exception as exc:                     # noqa: BLE001
        return (scan, acq, "", "", "", "", f"metadata read failed: {exc}")

    short, full = describe(scan)
    md["/entry/sample/short_description"] = short
    md["/entry/sample/full_description"] = full

    os.makedirs(out_dir, exist_ok=True)
    # Relink rather than copy: the frames are 27 MB each.
    if os.path.islink(out_bin) or os.path.exists(out_bin):
        os.remove(out_bin)
    os.symlink(src_bin, out_bin)

    # xpcs_schema.copy() is shallow and update_schema_at_runtime mutates in
    # place, so each acquisition needs its own deep copy.
    schema = update_schema_at_runtime(copy.deepcopy(xpcs_schema), md)
    create_nexus_format_metadata(out_md, schema)

    return (scan, acq, md["/entry/start_time"],
            f"{md['/entry/instrument/incident_beam/incident_beam_intensity']:.0f}",
            f"{md['/entry/sample/qnw3_temperature_set']:.2f}",
            " ".join(missing), "")


def main():
    # convert() reads these at module scope so the pool workers inherit them.
    global DEST, SCAN_GLOB

    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--scan-glob", default=SCAN_GLOB,
                    help="which scan groups to convert [default: %(default)s]")
    ap.add_argument("--dest", default=DEST,
                    help="output tree [default: %(default)s]")
    ap.add_argument("--manifest", default=MANIFEST,
                    help="where to record what was converted")
    args = ap.parse_args()
    DEST, SCAN_GLOB = args.dest, args.scan_glob

    jobs = []
    for scan_dir in sorted(glob.glob(os.path.join(ARCHIVE, args.scan_glob))):
        if not os.path.isdir(scan_dir):
            continue
        scan = os.path.basename(scan_dir)
        subs = [d for d in sorted(glob.glob(os.path.join(scan_dir, f"{scan}_*")))
                if os.path.isdir(d)]
        if subs:
            jobs.extend((scan, os.path.basename(d), d) for d in subs)
        elif os.path.exists(os.path.join(scan_dir, f"{scan}.bin")):
            # Some early scans are a single acquisition with the frame sitting
            # directly in the scan directory rather than in a sub-directory.
            jobs.append((scan.rsplit("_", 1)[0], scan, scan_dir))
    if args.limit:
        jobs = jobs[: args.limit]

    by_scan = {}
    for scan, _, _ in jobs:
        by_scan[scan] = by_scan.get(scan, 0) + 1
    print(f"{len(jobs)} acquisitions in {len(by_scan)} scan groups:")
    for scan, n in sorted(by_scan.items()):
        print(f"  {scan:34s} {n:5d}")

    if not args.apply:
        print(f"\ndry run; would write under {args.dest}")
        return

    os.makedirs(args.dest, exist_ok=True)
    with multiprocessing.Pool(args.procs) as pool:
        rows = pool.map(convert, jobs, chunksize=8)

    failures = [r for r in rows if r[6]]
    with open(args.manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["scan", "acquisition", "start_time", "I0Monitor", "qnw3_set",
                    "fields_not_recorded", "error"])
        w.writerows(rows)

    print(f"\nconverted {len(rows) - len(failures)} of {len(rows)}")
    print(f"manifest {args.manifest}")
    if failures:
        print(f"{len(failures)} failures, first few:")
        for r in failures[:10]:
            print("   ", r[1], r[6])
    partial = [r for r in rows if r[5] and not r[6]]
    if partial:
        print(f"{len(partial)} converted with fields the archive never recorded, e.g. {partial[0][5]}")


if __name__ == "__main__":
    main()

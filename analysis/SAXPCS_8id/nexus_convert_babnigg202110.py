"""Convert the 2021-3 isothermal series from Data Exchange to NeXus.

BEAMLINE-ONLY: needs /gdata (amber, adamite, amazonite) and the nexus_xpcs_aps
package.  Kept for provenance: it records how the NeXus tree was made.

For every acquisition of the H04* scan groups of babnigg202110 it writes

    babnigg202110_nexus/<acq>/<acq>.bin             symlink to the archive frames
    babnigg202110_nexus/<acq>/<acq>_metadata.hdf    NeXus metadata

and records each acquisition in a manifest.  The archive is opened read-only.
The ion chambers, ring current, start and end times, geometry, detector name
and sample temperatures are read from the archive; a field an acquisition did
not record is written as NaN, never as a placeholder.  The beam centre is the
refined (x=510.5, y=259.5) in NeXus orientation, and beam_center_position_x/y
equal position_x/y so the reader adds no offset; the README explains both.

    python nexus_convert_babnigg202110.py                 # dry run
    python nexus_convert_babnigg202110.py --apply
    python nexus_convert_babnigg202110.py --apply --limit 20     # smoke test
"""

import argparse
import copy
import csv
import glob
import multiprocessing
import os
import re
import uuid
from datetime import datetime

import h5py
import numpy as np

from nexus_xpcs_aps.nexus_schema import xpcs_schema
from nexus_xpcs_aps.nexus_utils import update_schema_at_runtime, create_nexus_format_metadata

ARCHIVE = "/gdata/s8id-dmdtn/2021-3/babnigg202110"
DEST = "/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus"
SCAN_GLOB = "H04*"
MANIFEST = os.path.join(os.path.dirname(os.path.abspath(__file__)),
                        "nexus_convert_babnigg202110_manifest.csv")

DX_TIME_FORMAT = "%a %b %d %H:%M:%S %Y"      # 'Sun Oct 17 15:38:32 2021'
OUT_TIME_FORMAT = "%Y-%m-%d %H:%M:%S"        # what nexus_read.py parses

# Geometry in NeXus orientation, refined in pySimpleMask against the measured scattering.
BEAM_CENTER_X, BEAM_CENTER_Y = 510.5, 259.5
DETECTOR_NAME = "rigaku500k"
CYCLE, PROPOSAL = "2021-3", "babnigg202110"

# Sample label and nominal temperature, read off the scan name.
SCAN_RE = re.compile(r"^(H\d{4})_(D\d+)_(\d{3})C10p_att00_Rq0$")

# The groups read_metadata() needs in a Data Exchange record.
REQUIRED_GROUPS = ("measurement/instrument/acquisition", "measurement/instrument/detector",
                   "measurement/instrument/source_begin", "measurement/instrument/source_end",
                   "measurement/sample")


def scalar(obj):
    return float(np.asarray(obj[()]).ravel()[0])


def text(obj):
    v = np.asarray(obj[()]).ravel()[0]
    if isinstance(v, bytes):
        return v.decode("utf-8")
    return str(v)


def optional(group, name, missing):
    """A measured scalar, or NaN if this acquisition did not record it (the
    Oct 9 scans have no source_end/TransmissionMonitor).  The names of the
    missing fields are appended to missing."""
    if name in group:
        return scalar(group[name])
    missing.append(name)
    return float("nan")


def dx_time(obj):
    """A Data Exchange time stamp rewritten in the NeXus format."""
    return datetime.strptime(text(obj).strip(), DX_TIME_FORMAT).strftime(OUT_TIME_FORMAT)


def dx_file(acq_dir, acq):
    """The Data Exchange metadata file of one acquisition, or None."""
    direct = os.path.join(acq_dir, f"{acq}_0001-100000.hdf")
    if os.path.exists(direct):
        return direct
    found = sorted(glob.glob(os.path.join(acq_dir, "*_0001-*.hdf")))
    if found:
        return found[0]
    return None


def record_defect(path):
    """Why the Data Exchange record at path cannot be converted, or '' if it can.

    A few acquisitions of the archive have incomplete records: no detector
    group, or a start or end time of 'None'.  They are listed in the manifest
    and not converted.
    """
    with h5py.File(path, "r") as hf:
        for group in REQUIRED_GROUPS:
            if group not in hf:
                return f"no {group} in the Data Exchange record"
        for which in ("source_begin", "source_end"):
            if text(hf[f"measurement/instrument/{which}/datetime"]).strip() == "None":
                return f"no {which} time in the Data Exchange record"
    return ""


def read_metadata(path, acq):
    """Map one acquisition's Data Exchange record onto NeXus paths.

    Returns (metadata, missing), where missing names the measured fields this
    acquisition did not record.
    """
    missing = []
    with h5py.File(path, "r") as hf:
        acqg = hf["measurement/instrument/acquisition"]
        det = hf["measurement/instrument/detector"]
        beg = hf["measurement/instrument/source_begin"]
        end = hf["measurement/instrument/source_end"]
        smp = hf["measurement/sample"]

        translation = np.asarray(smp["translation"][()]).ravel()
        md = {
            "/entry/entry_identifier": acq,
            "/entry/entry_identifier_uuid": uuid.uuid4().hex,
            "/entry/beamline": "APS-8-ID-I",
            "/entry/scan_number": int(scalar(acqg["specscan_data_number"])),
            "/entry/start_time": dx_time(beg["datetime"]),
            "/entry/end_time": dx_time(end["datetime"]),
            "/entry/user/cycle": CYCLE,
            "/entry/user/proposal_id": PROPOSAL,
            # not recorded in the archive; "unknown" rather than the schema's
            # placeholder name, which would read as real
            "/entry/user/name": "unknown",
            "/entry/user/email": "unknown",
            "/entry/user/institution": "unknown",

            "/entry/instrument/detector_1/detector_name": DETECTOR_NAME,
            "/entry/instrument/detector_1/count_time": scalar(det["exposure_time"]),
            "/entry/instrument/detector_1/frame_time": scalar(det["exposure_period"]),
            # Data Exchange records these in mm; NeXus wants m
            "/entry/instrument/detector_1/distance": scalar(det["distance"]) * 1e-3,
            "/entry/instrument/detector_1/x_pixel_size": scalar(det["x_pixel_size"]) * 1e-3,
            "/entry/instrument/detector_1/y_pixel_size": scalar(det["y_pixel_size"]) * 1e-3,
            "/entry/instrument/detector_1/beam_center_x": BEAM_CENTER_X,
            "/entry/instrument/detector_1/beam_center_y": BEAM_CENTER_Y,
            # equal positions, so the reader adds no offset to the beam centre
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
    header, sample, temp = m.groups()
    # The three digits are the zone-3 setpoint times ten, truncated: 282 is
    # 28.25 C, 287 is 28.75 C.  The readback in the file is authoritative.
    return sample, f"ELP sample {sample}, nominal hold {int(temp) / 10.0:.2f} C"


def convert(job):
    """Convert one acquisition.  Returns a manifest row; its last entry says
    why the acquisition was not converted, or is empty."""
    scan, acq, acq_dir, dest = job
    src_bin = os.path.join(acq_dir, f"{acq}.bin")
    out_dir = os.path.join(dest, acq)
    out_bin = os.path.join(out_dir, f"{acq}.bin")
    out_md = os.path.join(out_dir, f"{acq}_metadata.hdf")

    dx = dx_file(acq_dir, acq)
    if dx is None:
        return (scan, acq, "", "", "", "", "no Data Exchange metadata")
    if not os.path.exists(src_bin):
        return (scan, acq, "", "", "", "", "no .bin in the archive")
    defect = record_defect(dx)
    if defect:
        return (scan, acq, "", "", "", "", defect)

    md, missing = read_metadata(dx, acq)
    short, full = describe(scan)
    md["/entry/sample/short_description"] = short
    md["/entry/sample/full_description"] = full

    os.makedirs(out_dir, exist_ok=True)
    # link rather than copy: the frames are 27 MB each
    if os.path.islink(out_bin) or os.path.exists(out_bin):
        os.remove(out_bin)
    os.symlink(src_bin, out_bin)

    # update_schema_at_runtime changes the schema in place, so each
    # acquisition gets its own deep copy
    schema = update_schema_at_runtime(copy.deepcopy(xpcs_schema), md)
    create_nexus_format_metadata(out_md, schema)

    return (scan, acq, md["/entry/start_time"],
            f"{md['/entry/instrument/incident_beam/incident_beam_intensity']:.0f}",
            f"{md['/entry/sample/qnw3_temperature_set']:.2f}",
            " ".join(missing), "")


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--apply", action="store_true")
    ap.add_argument("--limit", type=int, default=0)
    ap.add_argument("--procs", type=int, default=16)
    ap.add_argument("--scan-glob", nargs="+", default=[SCAN_GLOB],
                    help="one or more scan-group patterns [default: %(default)s]")
    ap.add_argument("--dest", default=DEST,
                    help="output tree [default: %(default)s]")
    ap.add_argument("--manifest", default=MANIFEST,
                    help="where to record what was converted")
    args = ap.parse_args()

    scan_dirs = set()
    for pattern in args.scan_glob:
        for d in glob.glob(os.path.join(ARCHIVE, pattern)):
            scan_dirs.add(d)
    jobs = []
    for scan_dir in sorted(scan_dirs):
        if not os.path.isdir(scan_dir):
            continue
        scan = os.path.basename(scan_dir)
        acquisition_dirs = []
        for d in sorted(glob.glob(os.path.join(scan_dir, f"{scan}_*"))):
            if os.path.isdir(d):
                acquisition_dirs.append(d)
        if acquisition_dirs:
            for d in acquisition_dirs:
                jobs.append((scan, os.path.basename(d), d, args.dest))
        elif os.path.exists(os.path.join(scan_dir, f"{scan}.bin")):
            # a few early scans are one acquisition, the frames directly in the scan directory
            jobs.append((scan.rsplit("_", 1)[0], scan, scan_dir, args.dest))
    if args.limit:
        jobs = jobs[: args.limit]

    by_scan = {}
    for job in jobs:
        by_scan[job[0]] = by_scan.get(job[0], 0) + 1
    print(f"{len(jobs)} acquisitions in {len(by_scan)} scan groups:")
    for scan, n in sorted(by_scan.items()):
        print(f"  {scan:34s} {n:5d}")

    if not args.apply:
        print(f"\ndry run; would write under {args.dest}")
        return

    os.makedirs(args.dest, exist_ok=True)
    with multiprocessing.Pool(args.procs) as pool:
        rows = pool.map(convert, jobs, chunksize=8)

    with open(args.manifest, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["scan", "acquisition", "start_time", "I0Monitor", "qnw3_set",
                    "fields_not_recorded", "error"])
        w.writerows(rows)

    failures = []
    partial = []
    for r in rows:
        if r[6]:
            failures.append(r)
        elif r[5]:
            partial.append(r)
    print(f"\nconverted {len(rows) - len(failures)} of {len(rows)}")
    print(f"manifest {args.manifest}")
    if failures:
        print(f"{len(failures)} failures, first few:")
        for r in failures[:10]:
            print("   ", r[1], r[6])
    if partial:
        print(f"{len(partial)} converted with fields the archive never recorded, e.g. {partial[0][5]}")


if __name__ == "__main__":
    main()

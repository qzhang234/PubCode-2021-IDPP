#!/bin/bash
# Correlate every acquisition of the 2021-3 isothermal series.
#
# BEAMLINE-ONLY. Needs /gdata and a GPU host (adamite, amazonite). Written to
# run during the 5-6 Oct 2026 APS shutdown, and kept here for provenance.
#
# Scheduled with `at` on the GPU host rather than from a Claude session, so the
# run does not depend on anything staying alive:
#     echo /path/to/run_babnigg202110_correlations.sh | at 08:00 2026-10-05
#
# Design notes
# ------------
# -i -2 takes ONE GPU per process, by file lock in /dev/shm/gpu_locks; it does
# not spread a single invocation over the card set. Parallelism therefore comes
# from running NWORKERS invocations at once, which the scheduler spreads over
# the free GPUs and queues when they are busy. Measured 0.70 files/s at 4
# workers on adamite's four A100s, and 8 workers gave 0.707, so the GPUs are
# saturated at four and there is nothing to gain by oversubscribing.
#
# --crop-ratio-threshold defaults to 1.0, which crops whenever any pixel is
# masked, and cropping makes every g2 come back NaN on this data. 0.5 keeps
# cropping off for this qmap, whose valid fraction is 0.87.
#
# -w overwrites in place. Without it boost_corr appends _1 to the name of any
# result that already exists, which is how duplicates appeared before.
#
# DEADLINE is a safety net for the beam coming back, not the expected finish.
# No new chunk starts after it; chunks already running are left to finish,
# which takes under a minute. The whole run is about three hours.

set -u

HOST=$(hostname -s)
NEXUS=/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus
QMAP=$NEXUS/babnigg202110_nexus_Sq270_Dq27_Sphi8_Dphi1_Q0p0032_Lin_CrossRm.hdf
OUT=$NEXUS/reprocess_results
RUNDIR=$NEXUS/correlation_run_20261005
NWORKERS=${NWORKERS:-4}
CHUNK=${CHUNK:-12}
# 2026-10-06 07:00 local, an hour before the beam returns.
DEADLINE=${DEADLINE:-1791288000}

export PATH=/home/beams/8IDIUSER/bin:$PATH
mkdir -p "$RUNDIR" "$OUT"
LOG=$RUNDIR/progress_${HOST}.log
STATE=$RUNDIR/state_${HOST}

log() { echo "$(date '+%F %T') $*" | tee -a "$LOG"; }

log "=== start on $HOST, workers=$NWORKERS chunk=$CHUNK deadline=$(date -d @$DEADLINE '+%F %T')"
[ -f "$QMAP" ] || { log "FATAL: qmap missing: $QMAP"; exit 1; }
log "qmap $(basename "$QMAP")"

# Any _1 files are leftovers from a run that lacked -w; they would shadow the
# real results in a listing, so clear them before starting.
nd=$(ls "$OUT"/*_results_1.hdf 2>/dev/null | wc -l)
[ "$nd" -gt 0 ] && { rm -f "$OUT"/*_results_1.hdf; log "removed $nd stale _1 files"; }

LIST=$RUNDIR/filelist.txt
if [ ! -s "$LIST" ]; then
    find "$NEXUS" -mindepth 1 -maxdepth 1 -type d ! -name 'reprocess_results' \
         ! -name 'correlation_run_*' -printf '%f\n' | sort | while read -r a; do
        [ -e "$NEXUS/$a/$a.bin" ] && echo "$NEXUS/$a/$a.bin"
    done > "$LIST"
fi
TOTAL=$(wc -l < "$LIST")
log "acquisitions to process: $TOTAL"

CHUNKDIR=$RUNDIR/chunks_${HOST}
rm -rf "$CHUNKDIR"; mkdir -p "$CHUNKDIR"
split -l "$CHUNK" -d -a 4 "$LIST" "$CHUNKDIR/c_"
NCHUNK=$(ls "$CHUNKDIR" | wc -l)
log "split into $NCHUNK chunks of $CHUNK"

run_chunk() {
    c="$1"
    if [ "$(date +%s)" -ge "$DEADLINE" ]; then return 0; fi
    boost_corr_bin -r $(tr '\n' ' ' < "$c") -q "$QMAP" -o "$OUT" \
        -i -2 -t Multitau --crop-ratio-threshold 0.5 -w \
        >> "$RUNDIR/boost_${HOST}.log" 2>&1
    echo "$(date '+%F %T') $(basename "$c") rc=$?" >> "$STATE"
}
export -f run_chunk
export QMAP OUT RUNDIR HOST STATE DEADLINE PATH

T0=$(date +%s)
ls "$CHUNKDIR"/c_* | xargs -r -P "$NWORKERS" -I{} bash -c 'run_chunk "$@"' _ {}
T1=$(date +%s)

DONE=$(ls "$OUT"/*_results.hdf 2>/dev/null | wc -l)
DUPE=$(ls "$OUT"/*_results_1.hdf 2>/dev/null | wc -l)
log "=== finished in $(( (T1-T0)/60 )) min: $DONE results present, $DUPE duplicates"

# Name anything that did not produce a result, so a short follow-up run can
# pick it up rather than repeating the whole set.
MISSING=$RUNDIR/missing_${HOST}.txt
: > "$MISSING"
while read -r b; do
    a=$(basename "$b" .bin)
    [ -e "$OUT/${a}_results.hdf" ] || echo "$b" >> "$MISSING"
done < "$LIST"
log "missing: $(wc -l < "$MISSING") (listed in $(basename "$MISSING"))"

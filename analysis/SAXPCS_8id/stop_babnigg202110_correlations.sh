#!/bin/bash
# Hard-stop the 2021-3 correlation run before the beam returns.
#
# BEAMLINE-ONLY. Scheduled with `at` on the GPU host for 07:30 on 6 Oct 2026,
# an hour before the beam comes back. The run itself is expected to finish
# around 11:00 on 5 Oct, so this should normally find nothing to do.
#
# The pattern lives in this file rather than in the `at` command because
# pkill -f matches command lines: an `at` job that carried the pattern inline
# would match its own shell and kill itself before logging. A script's command
# line is just its path, so there is nothing for the pattern to match.
#
# Only processes whose command line names this cycle's NeXus tree are killed,
# so unrelated boost_corr work by other groups on the same host is left alone.

set -u
PATTERN='boost_corr.*babnigg202110_nexus'
RUNDIR=/gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/correlation_run_20261005
LOG=$RUNDIR/progress_$(hostname -s).log
mkdir -p "$RUNDIR"

# Never signal this script, its shell, or its ancestors.
SELF=$$
EXCLUDE=" $SELF $PPID "
p=$PPID
while [ -n "$p" ] && [ "$p" -gt 1 ] 2>/dev/null; do
    p=$(ps -o ppid= -p "$p" 2>/dev/null | tr -d ' ')
    [ -n "$p" ] && EXCLUDE="$EXCLUDE $p "
done

targets=$(pgrep -u "$(id -u)" -f "$PATTERN" 2>/dev/null | while read -r pid; do
    case "$EXCLUDE" in *" $pid "*) continue ;; esac
    echo "$pid"
done)

if [ -z "$targets" ]; then
    echo "$(date '+%F %T') hard stop: nothing running, run had already finished" >> "$LOG"
    exit 0
fi

echo "$(date '+%F %T') hard stop before beam return, signalling: $(echo $targets | tr '\n' ' ')" >> "$LOG"
kill -TERM $targets 2>/dev/null
sleep 20
still=$(for pid in $targets; do kill -0 "$pid" 2>/dev/null && echo "$pid"; done)
[ -n "$still" ] && { kill -KILL $still 2>/dev/null; echo "$(date '+%F %T') SIGKILL sent to $(echo $still | tr '\n' ' ')" >> "$LOG"; }

printf 'To: %s\nFrom: 8idiuser@%s.xray.aps.anl.gov\nSubject: %s\n\n%s\n' \
    "qzhang234@anl.gov" "$(hostname -s)" \
    "babnigg202110 correlation run hard-stopped on $(hostname -s) before beam return" \
"The 07:30 safety stop found the run still going and ended it so the GPUs are
free for the beam at 08:00.

results present : $(ls /gdata/s8id-dmdtn/2021-3/babnigg202110_nexus/reprocess_results/*_results.hdf 2>/dev/null | wc -l) of 7421
run directory   : $RUNDIR

$(tail -15 "$LOG" 2>/dev/null)" | /usr/sbin/sendmail -t -i 2>/dev/null

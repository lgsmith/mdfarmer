#!/usr/bin/env bash
# Start, stop and check on the shakedown-omm tender.
#
# The tender is the long-lived python process that submits this campaign's
# Slurm jobs and re-submits them as generations finish. It runs detached on the
# machine you launch it from -- a login node or a workstation -- and NOT inside
# an sbatch job of its own: it sleeps between ticks, so a whole allocation would
# sit idle holding it, and that allocation's walltime (or its preemption) would
# end the campaign with it. All it needs from the cluster is sbatch.
#
#   ./drive_omm.sh --check         # run shape and readiness; writes nothing
#   ./drive_omm.sh --dry-run       # dirs, configs, job scripts; submits nothing
#   ./drive_omm.sh                 # start the tender, detached, and return
#   ./drive_omm.sh --status        # up or down, its pid, the tail of its log
#   ./drive_omm.sh --stop          # brake it at its next tick
#
# Arguments after an explicit `start`, `--check` or `--dry-run` are passed on to
# farmer.py: `./drive_omm.sh start --n-gens 1`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
PROJECT='shakedown-omm'
CAMPAIGN="$HERE/data/$PROJECT"       # farmer.py's TRAJ_TOP
LOG="$HERE/$PROJECT.tend.out"
LOCK="$CAMPAIGN/tender.lock"         # one tender per campaign, held while it runs
STOP="$HERE/stop"                    # Farmer brakes on ./stop in its cwd

# `mamba run --no-capture-output` is broken here, so python -u is what keeps
# the log live.
CONDA_ENV="${CONDA_ENV:-omm}"
RUN_PY=(mamba run -n "$CONDA_ENV" python -u)

GAP="${GAP:-60}"                     # seconds before re-entering a failed tender
MIN_RUN_S=60                         # a shorter run is a broken setup, not a hiccup
MAX_FAST_EXITS=3
START_TIMEOUT_S=30                   # how long to wait for the detached tender

cd "$HERE"


usage() {
    cat <<EOF
usage: $(basename "$SELF") [--check|--dry-run|start|--status|--stop] [farmer.py args]

  --check     run shape and input readiness; writes nothing, starts nothing
  --dry-run   every directory, config and job script; submits nothing
  start       start the tender detached (the default with no argument)
  --status    whether a tender is up, its pid, and the tail of its log
  --stop      write the brake file; the tender stops at its next tick
EOF
}


tender_running() {
    mkdir -p "$(dirname "$LOCK")"
    ! flock -n "$LOCK" true 2>/dev/null
}


# What the lock file says; stale unless tender_running agrees.
last_tender() {
    if [ -s "$LOCK" ]; then cat "$LOCK"; else echo "none recorded"; fi
}


preflight() {
    if ! command -v mamba >/dev/null; then
        echo "no mamba on PATH; cannot reach the $CONDA_ENV environment" >&2
        exit 1
    fi
    if ! command -v sbatch >/dev/null; then
        echo "no sbatch on PATH; run this where jobs can be submitted" >&2
        exit 1
    fi
    # Which checkout the campaign will actually run, since `import mdfarmer`
    # resolves to whatever the environment installed, not necessarily this tree.
    "${RUN_PY[@]}" -c 'import mdfarmer; print("mdfarmer:", mdfarmer.__file__)'
}


plan() {
    echo "project   $PROJECT"
    echo "cwd       $HERE"
    echo "tender    ${RUN_PY[*]} farmer.py $*"
    echo "log       $LOG (appended)"
    echo "lock      $LOCK"
    echo "brake     $STOP"
    echo "re-enter  every ${GAP}s until the campaign finishes or is braked"
}


start_tender() {
    if tender_running; then
        echo "REFUSING: a tender for $PROJECT already holds $LOCK" >&2
        echo "  it says: $(last_tender)" >&2
        echo "  log: $LOG" >&2
        echo "  stop it with: $SELF --stop" >&2
        exit 1
    fi
    preflight
    if [ -e "$STOP" ]; then
        echo "clearing the leftover brake file $STOP"
        rm -f "$STOP"
    fi
    plan "$@"
    mkdir -p "$(dirname "$LOCK")"
    setsid nohup bash "$SELF" --loop "$@" >>"$LOG" 2>&1 </dev/null &
    local waited=0
    while [ "$waited" -lt "$START_TIMEOUT_S" ]; do
        if tender_running; then
            echo "tender up: $(last_tender)"
            echo "follow it with: tail -f $LOG"
            return 0
        fi
        sleep 1
        waited=$((waited + 1))
    done
    echo "the tender did not take $LOCK within ${START_TIMEOUT_S}s; last log lines:" >&2
    tail -n 20 "$LOG" >&2 || true
    exit 1
}


# The re-entering loop, run detached by start_tender. Farmer.launch drops a
# clone for good on one transient sbatch failure, and a fresh tender rebuilds
# every clone from disk and re-adopts the running job ids, so re-entering is
# the recovery.
tender_loop() {
    exec 9>>"$LOCK"
    if ! flock -n 9; then
        echo "another tender holds $LOCK; this one exits"
        exit 1
    fi
    printf 'pid %s on %s since %s\n' "$$" "$(hostname)" "$(date -Is)" >"$LOCK"
    local fast_exits=0 started status
    while [ ! -e "$STOP" ]; do
        echo "=== tender starting $(date -Is): ${RUN_PY[*]} farmer.py $* ==="
        started=$SECONDS
        status=0
        "${RUN_PY[@]}" farmer.py "$@" || status=$?
        echo "=== tender exited $status at $(date -Is) ==="
        if [ "$status" -eq 0 ]; then
            echo "=== every clone finished; not re-entering ==="
            break
        fi
        if [ -e "$STOP" ]; then
            break
        fi
        if [ $((SECONDS - started)) -lt "$MIN_RUN_S" ]; then
            fast_exits=$((fast_exits + 1))
        else
            fast_exits=0
        fi
        if [ "$fast_exits" -ge "$MAX_FAST_EXITS" ]; then
            echo "=== $fast_exits tenders in a row died inside ${MIN_RUN_S}s;" \
                 "that is a setup error, not a scheduler hiccup. Giving up. ==="
            exit 1
        fi
        echo "=== re-entering in ${GAP}s ==="
        sleep "$GAP"
    done
    if [ -e "$STOP" ]; then
        echo "=== brake file $STOP is present; loop done $(date -Is) ==="
    fi
    echo "=== jobs already submitted keep running; squeue -u \$USER to see them ==="
}


show_status() {
    if tender_running; then
        echo "tender UP    $(last_tender)"
    else
        echo "tender DOWN  (last: $(last_tender))"
    fi
    echo "log          $LOG"
    echo "brake        $STOP$([ -e "$STOP" ] && echo ' (present)')"
    if [ -f "$LOG" ]; then
        echo "--- last 5 log lines ---"
        tail -n 5 "$LOG"
    fi
}


stop_tender() {
    if ! tender_running; then
        echo "no tender holds $LOCK; nothing to brake"
        exit 0
    fi
    touch "$STOP"
    echo "wrote the brake file $STOP"
    echo "the tender stops at its next tick, and the loop then exits;"
    echo "jobs already submitted keep running. Watch: tail -f $LOG"
}


mode="${1:---start}"
if [ "$#" -gt 0 ]; then shift; fi
case "$mode" in
    --start|start) start_tender "$@" ;;
    --loop)        tender_loop "$@" ;;
    --check)       preflight; "${RUN_PY[@]}" farmer.py --check "$@" ;;
    --dry-run)     preflight; plan --dry-run "$@"
                   "${RUN_PY[@]}" farmer.py --dry-run "$@" ;;
    --status)      show_status ;;
    --stop)        stop_tender ;;
    -h|--help)     usage ;;
    *)             usage >&2; exit 2 ;;
esac

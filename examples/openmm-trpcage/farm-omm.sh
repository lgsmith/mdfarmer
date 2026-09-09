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
#   ./farm-omm.sh --check         # run shape and readiness; writes nothing
#   ./farm-omm.sh --dry-run       # dirs, configs, job scripts; submits nothing
#   ./farm-omm.sh                 # start the tender, detached, and return
#   ./farm-omm.sh --status        # up or down, its pid, the tail of its log
#   ./farm-omm.sh --stop          # brake it at its next tick
#
# The conda environment holding mdfarmer comes from --env, else $CONDA_ENV,
# else a .conda-env file at the top of the checkout.
# Arguments after an explicit `start`, `--check` or `--dry-run` are passed on to
# farmer.py: `./farm-omm.sh start --n-gens 1`.
set -euo pipefail

HERE="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SELF="$HERE/$(basename "${BASH_SOURCE[0]}")"
PROJECT='shakedown-omm'
CAMPAIGN="$HERE/data/$PROJECT"       # farmer.py's TRAJ_TOP
LOG="$HERE/$PROJECT.tend.out"
LOCK="$CAMPAIGN/tender.lock"         # one tender per campaign, held while it runs
STOP="$HERE/stop"                    # Farmer brakes on ./stop in its cwd
BRAKED_EXIT=2                        # farmer.py's status for "asked to stop"

# The environment holding mdfarmer and its engine. Named, not hard-coded: this
# example is meant to be run by someone whose env is not called what mine is.
# Resolution order: --env, then $CONDA_ENV, then a .conda-env file at the top of
# the checkout. With one of those in place, launching is just ./$(basename "$0").
CONDA_ENV_FILE="$HERE/../../.conda-env"
CONDA_ENV="${CONDA_ENV:-}"
RUN_PY=(python -u)               # what activate_env leaves on PATH

GAP="${GAP:-60}"                     # seconds before re-entering a failed tender
MIN_RUN_S=60                         # a shorter run is a broken setup, not a hiccup
MAX_FAST_EXITS=3
START_TIMEOUT_S=30                   # how long to wait for the detached tender

cd "$HERE"


usage() {
    cat <<EOF
usage: $(basename "$SELF") [--env NAME] [--check|--dry-run|start|--status|--stop] [farmer.py args]

  --env NAME  conda env holding mdfarmer; else \$CONDA_ENV, else .conda-env
  --check     run shape and input readiness; writes nothing, starts nothing
  --dry-run   every directory, config and job script; submits nothing
  start       start the tender detached (the default with no argument)
  --status    whether a tender is up, its pid, and the tail of its log
  --stop      write the brake file; the tender stops at its next tick
EOF
}


# Put CONDA_ENV's bin on PATH, so `python` below is a direct child of the loop:
# its exit status is the one the loop reads and its stdout is not buffered by a
# wrapper. set +u because conda's own shell functions do not survive it.
activate_env() {
    if [ -z "$CONDA_ENV" ] && [ -r "$CONDA_ENV_FILE" ]; then
        CONDA_ENV="$(tr -d '[:space:]' < "$CONDA_ENV_FILE")"
    fi
    if [ -z "$CONDA_ENV" ]; then
        echo "no environment named. Pass one:" >&2
        echo "  $SELF --env NAME [--check|--dry-run|start|...]" >&2
        echo "or set CONDA_ENV=NAME, or write the name into" >&2
        echo "  $CONDA_ENV_FILE" >&2
        echo "It must hold mdfarmer and its engine." >&2
        exit 2
    fi
    local base
    base="$(conda info --base 2>/dev/null)" || {
        echo "no conda on PATH; cannot activate $CONDA_ENV" >&2; exit 1; }
    set +u
    # shellcheck disable=SC1091
    source "$base/etc/profile.d/conda.sh"
    conda activate "$CONDA_ENV" || {
        echo "could not activate the environment '$CONDA_ENV'" >&2; exit 1; }
    set -u
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
    if ! command -v sbatch >/dev/null; then
        echo "no sbatch on PATH; run this where jobs can be submitted" >&2
        exit 1
    fi
    activate_env
    echo "env       $CONDA_ENV ($(command -v python))"
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
    mkdir -p "$(dirname "$LOCK")"
    exec 9>>"$LOCK"
    if ! flock -n 9; then
        echo "another tender holds $LOCK; this one exits"
        exit 1
    fi
    activate_env
    printf 'pid %s on %s since %s\n' "$$" "$(hostname)" "$(date -Is)" >"$LOCK"
    local fast_exits=0 started status
    while true; do
        echo "=== tender starting $(date -Is): ${RUN_PY[*]} farmer.py $* ==="
        started=$SECONDS
        status=0
        "${RUN_PY[@]}" farmer.py "$@" || status=$?
        echo "=== tender exited $status at $(date -Is) ==="
        if [ "$status" -eq 0 ]; then
            echo "=== every clone finished; not re-entering ==="
            break
        fi
        # The tender distinguishes "asked to stop" from "died", so this loop
        # reads its status and forms no opinion about the brake file itself.
        if [ "$status" -eq "$BRAKED_EXIT" ]; then
            echo "=== tender was braked; not re-entering ==="
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


# --env NAME comes before the mode, and is exported so the detached loop and
# every re-entered tender inherit it.
while [ "$#" -gt 0 ]; do
    case "$1" in
        --env)   CONDA_ENV="${2:-}"; shift 2 ;;
        --env=*) CONDA_ENV="${1#--env=}"; shift ;;
        *)       break ;;
    esac
done
export CONDA_ENV

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

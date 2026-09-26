#!/usr/bin/env bash
# Single entry point.  ./run.sh <run_id> [nextflow args...]
# Pins the launch dir, workDir, outdir and log so a run can never scatter.
set -euo pipefail

RUN_ID="${1:?usage: ./run.sh <run_id> [nextflow args...]}"; shift
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
RUN_DIR="$PIPE_DIR/../Runs/$RUN_ID"
# work/ is disposable and I/O-heavy -> local NVMe, not the CephFS home.
SCRATCH="${NXF_SCRATCH_ROOT:-/tmp/$USER-nf}"
WORK_DIR="$SCRATCH/$RUN_ID"

# A mistyped --param is silently null in Nextflow. Reject anything no config declares.
DECLARED=$(grep -rhoE '^[[:space:]]*[a-zA-Z_][a-zA-Z0-9_]*[[:space:]]*=' \
             "$PIPE_DIR/nextflow.config" "$PIPE_DIR"/conf/*.config | tr -d ' =' | sort -u)
for a in "$@"; do
  case "$a" in
    --*) n="${a#--}"
         grep -qx "$n" <<<"$DECLARED" || { echo "unknown parameter: --$n" >&2; exit 1; } ;;
  esac
done

# Rtmp must exist or R silently falls back to /tmp and the setting is lost.
mkdir -p "$RUN_DIR/results" "$RUN_DIR/logs" "$WORK_DIR" "$SCRATCH/Rtmp"
ln -sfn "$WORK_DIR" "$RUN_DIR/work"   # failed task dirs stay reachable from Runs/

# Machine-specific paths live in conf/site.yaml, outside git (see conf/site.example.yaml)
SITE=()
[ -f "$PIPE_DIR/conf/site.yaml" ] && SITE=(-params-file "$PIPE_DIR/conf/site.yaml")

printf 'run_id  : %s\nresults : %s\nwork    : %s (local scratch, NOT backed up)\nlog     : %s\nsite    : %s\n\n' \
  "$RUN_ID" "$RUN_DIR/results" "$WORK_DIR" "$RUN_DIR/logs/nextflow.log" \
  "$([ ${#SITE[@]} -gt 0 ] && echo conf/site.yaml || echo 'none (conf/site.yaml not found)')"

cd "$PIPE_DIR"
nextflow -log "$RUN_DIR/logs/nextflow.log" \
    run main.nf \
    ${SITE[@]+"${SITE[@]}"} \
    -work-dir "$WORK_DIR" \
    --outdir  "$RUN_DIR/results" \
    --run_id  "$RUN_ID" \
    "$@"

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

# A bare -resume picks the last run in this folder, whatever its stage, and each stage's
# cached tasks live in its own session. Resume this run's last session of the same --stage
# instead; with none, start fresh. An explicit session id is kept.
# Nextflow does not notice a stage's input object (a saveGiotto folder) being rewritten, so
# after an earlier stage reran, a resume would reuse tasks built on the old object. When
# nothing the stage wrote is newer than its input object, start fresh instead. Cached runs
# rewrite neither, so the history's run times cannot tell.
ARGS=("$@")
STAGE=""
BASIS=""
IN=""
for ((i = 0; i < ${#ARGS[@]}; i++)); do
  case "${ARGS[i]}" in
    --stage) STAGE="${ARGS[i+1]:-}" ;;
    --cluster_basis) BASIS="${ARGS[i+1]:-}" ;;
    --merged_giotto) IN="${ARGS[i+1]:-}" ;;
  esac
done
# Without --cluster_basis, annotate's basis comes from the run's labels folder, as main.nf
# reads it: cluster_basis.txt, or the one cluster_labels.<basis>.csv there
if [ "$STAGE" = annotate ] && [ -z "$BASIS" ]; then
  LBL="$RUN_DIR/labels"
  if [ -f "$LBL/cluster_basis.txt" ]; then
    BASIS=$(head -n 1 "$LBL/cluster_basis.txt" | tr -d '[:space:]')
  else
    FOUND=("$LBL"/cluster_labels.*.csv)
    if [ ${#FOUND[@]} -eq 1 ] && [ -f "${FOUND[0]}" ]; then
      BASIS=$(basename "${FOUND[0]}" .csv)
      BASIS=${BASIS#cluster_labels.}
    fi
  fi
fi
# The object each stage reads, as main.nf picks it, and where the stage writes. Stages that
# rewrite the object they read (subcluster, supervised) or read raw data are not checked
OUT=""
case "$STAGE" in
  annotate)
    OUT=objects/annotate/giotto
    if [[ "$BASIS" == sub_* ]]; then DEF=subcluster; else DEF=discover; fi ;;
  spatial) OUT=objects/spatial/giotto; DEF=annotate ;;
  cci) OUT=tables/08_cci; DEF=spatial ;;
  de) OUT=tables/09_de; DEF=spatial ;;
  banksy_sweep) OUT=tables/07_spatial/banksy_sweep; DEF=spatial ;;
esac
if [ -n "$OUT" ]; then
  if [ -z "$IN" ]; then
    IN="$RUN_DIR/results/objects/$DEF/giotto"
  elif [[ "$IN" != /* ]]; then
    IN="$PIPE_DIR/$IN"
  fi
fi
HIST="$PIPE_DIR/.nextflow/history"
# True when the input object exists and no file the stage wrote is newer than it
input_newer() {
  [ -n "$OUT" ] && [ -f "$IN/gobject.RDS" ] || return 1
  [ -z "$(find "$RUN_DIR/results/$OUT" -type f -newer "$IN/gobject.RDS" -print -quit 2>/dev/null || true)" ]
}
RESUME_NOTE=""
for ((i = 0; i < ${#ARGS[@]}; i++)); do
  [ "${ARGS[i]}" = "-resume" ] || continue
  NEXT="${ARGS[i+1]:-}"
  if [[ "$NEXT" =~ ^[0-9a-f]{8}(-[0-9a-f]{4}){3}-[0-9a-f]{12}$ ]]; then
    HIT=$(awk -F'\t' -v id="$NEXT" '$6 == id { nm = $3; t = $1 }
      END { if (t != "") print id "\t" nm "\t" t }' "$HIST" 2>/dev/null || true)
    [ -n "$HIT" ] || break
    TAKE=2
  else
    [ -n "$STAGE" ] || break
    HIT=$(awk -F'\t' -v r="$RUN_ID" -v s="$STAGE" '
      $NF ~ ("--run_id " r "( |$)") && $NF ~ ("--stage " s "( |$)") { id = $6; nm = $3; t = $1 }
      END { if (id != "") print id "\t" nm "\t" t }' "$HIST" 2>/dev/null || true)
    TAKE=1
    if [ -z "$HIT" ]; then
      ARGS=("${ARGS[@]:0:i}" "${ARGS[@]:i+1}")
      RESUME_NOTE="resume  : no earlier $STAGE run for $RUN_ID, starting fresh"
      break
    fi
  fi
  IFS=$'\t' read -r SID NAME _ <<<"$HIT"
  if input_newer; then
    ARGS=("${ARGS[@]:0:i}" "${ARGS[@]:i+TAKE}")
    RESUME_NOTE="resume  : ${IN#"$RUN_DIR/"} is newer than every $STAGE output, starting fresh"
  elif [ "$TAKE" = 1 ]; then
    ARGS=("${ARGS[@]:0:i+1}" "$SID" "${ARGS[@]:i+1}")
    RESUME_NOTE="resume  : last $STAGE run ($NAME), session $SID"
  else
    RESUME_NOTE="resume  : session $SID ($NAME)"
  fi
  break
done
set -- "${ARGS[@]+"${ARGS[@]}"}"

# Rtmp must exist or R silently falls back to /tmp and the setting is lost.
mkdir -p "$RUN_DIR/results" "$RUN_DIR/logs" "$WORK_DIR" "$SCRATCH/Rtmp"
ln -sfn "$WORK_DIR" "$RUN_DIR/work"   # failed task dirs stay reachable from Runs/

# Machine-specific paths live in conf/site.yaml, outside git (see conf/site.example.yaml)
SITE=()
[ -f "$PIPE_DIR/conf/site.yaml" ] && SITE=(-params-file "$PIPE_DIR/conf/site.yaml")

printf 'run_id  : %s\nresults : %s\nwork    : %s (local scratch, NOT backed up)\nlog     : %s\nsite    : %s\n\n' \
  "$RUN_ID" "$RUN_DIR/results" "$WORK_DIR" "$RUN_DIR/logs/nextflow.log" \
  "$([ ${#SITE[@]} -gt 0 ] && echo conf/site.yaml || echo 'none (conf/site.yaml not found)')"
[ -n "$RESUME_NOTE" ] && printf '%s\n\n' "$RESUME_NOTE"

cd "$PIPE_DIR"
nextflow -log "$RUN_DIR/logs/nextflow.log" \
    run main.nf \
    ${SITE[@]+"${SITE[@]}"} \
    -work-dir "$WORK_DIR" \
    --outdir  "$RUN_DIR/results" \
    --run_id  "$RUN_ID" \
    "$@"

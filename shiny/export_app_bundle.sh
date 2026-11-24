#!/usr/bin/env bash
# Exports a Giotto object to local_export/, the bundle that load_cosmx_nf.R reads, with bin/export_bundle.R
# The annotate stage writes the same bundle itself. This exports it from another object
#
#   ./export_app_bundle.sh [run_id] [genes] [poly_compress] [obj_subdir]
#   cosmx-nf <run_id> bundle [genes] [poly_compress] [obj_subdir]   (the same, in a project too)
# defaults: run_id=full  genes=all  poly_compress=xz  obj_subdir=objects/discover
# genes=N exports only the canonical markers and the top N DE genes per cluster
# poly_compress=none reports the vertex count and writes no polygons
# obj_subdir=objects/subcluster adds sub_<res> columns, objects/annotate also cell_type
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN=${1:-full}
GENES=${2:-all}
PCZ=${3:-xz}
OBJ=${4:-objects/discover}
# Runs live beside the pipeline, or in the project cosmx-nf names in COSMX_NF_RUNS
RUNS=${COSMX_NF_RUNS:-$PIPE/../Runs}
RES="$RUNS/$RUN/results"
[ -d "$RES/$OBJ/giotto" ] || { echo "no object at $RES/$OBJ/giotto" >&2; exit 1; }
LOG="$RUNS/$RUN/logs/export_app_bundle.txt"
mkdir -p "$(dirname "$LOG")" "$RES/local_export"

# The image and settings the stages use
"$PIPE/Environment/in_image.sh" Rscript --vanilla "$PIPE/bin/export_bundle.R" --input "$RES/$OBJ/giotto" \
  --outdir "$RES/local_export" --log "$LOG" --genes "$GENES" --poly_compress "$PCZ" \
  --markers "$PIPE/assets/canonical_markers.csv" --tables "$RES/tables"
echo "log: $LOG"

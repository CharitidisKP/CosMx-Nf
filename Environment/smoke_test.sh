#!/usr/bin/env bash
# Runs the smoke sheet with the current setup and with the image, then compares the tables
#   Environment/smoke_test.sh <image.sif> [assets] [smoke_old | smoke_new]
set -uo pipefail
# Ctrl-c stops the whole test, not only the current launch
trap 'echo "smoke_test.sh: stopped" >&2; exit 130' INT TERM

SIF=$(realpath "${1:?usage: Environment/smoke_test.sh <image.sif> [assets]}")
PIPE=$(cd "$(dirname "$0")/.." && pwd)
ASSETS=$(realpath "${2:-$PIPE/assets}")
RUNS=$(realpath -m "$PIPE/../Runs")
SCRATCH="${NXF_SCRATCH_ROOT:-/tmp/$USER-nf}"
BASIS=leiden_clus_res0.3
MAKE=(smoke_old smoke_new)
case ${3:-} in
  "") ;;
  smoke_old | smoke_new) MAKE=("$3") ;;
  *) echo "smoke_test.sh: run is smoke_old or smoke_new, not $3" >&2; exit 1 ;;
esac

[ -f "$SIF" ] || { echo "smoke_test.sh: no image at $SIF" >&2; exit 1; }
for d in profiles nichenet; do
  [ -d "$ASSETS/$d" ] || { echo "smoke_test.sh: no $d/ in $ASSETS" >&2; exit 1; }
done
# Each run made starts from an empty folder
for id in smoke_old smoke_new; do
  if [[ " ${MAKE[*]} " != *" $id "* ]]; then
    [ -d "$RUNS/$id/results" ] || { echo "smoke_test.sh: no $RUNS/$id to compare with" >&2; exit 1; }
  elif [ -e "$RUNS/$id" ] || [ -e "$SCRATCH/$id" ]; then
    echo "smoke_test.sh: $RUNS/$id or $SCRATCH/$id exists; move or delete it first" >&2
    exit 1
  fi
done

COMMON=(-profile server --input "$PIPE/assets/samplesheet.smoke.csv"
        --profile_dir "$ASSETS/profiles" --nichenet_model_dir "$ASSETS/nichenet")
# R with only the image's own library
IMAGE=(--sif "$SIF" --r_libs /none --r_libs_user /none
       --python_path /opt/py/bin/python --container_python /opt/py/bin/python)

LIBS=$(apptainer exec --cleanenv --env R_LIBS=/none --env R_LIBS_USER=/none "$SIF" \
  Rscript -e 'writeLines(.libPaths())') || exit 1
echo "R libraries the image's tasks will see:"
printf '  %s\n' $LIBS
case $LIBS in
  *"$HOME"*) echo "smoke_test.sh: a library under $HOME reaches the image's R" >&2; exit 1 ;;
esac

FAILED=()
# run <id> <stage> <what> [cosmx-nf args]: one fresh launch, noted in FAILED when it fails
run() {
  local id=$1 stage=$2 what=$3
  shift 3
  echo
  echo "=== $id: $what"
  "$PIPE/cosmx-nf" "$id" "$stage" --fresh "${COMMON[@]}" "$@" || { FAILED+=("$id: $what"); return 1; }
}

# Labels each cluster with the diagnostics' best guess
draft_labels() {
  local r="$RUNS/$1" f="cluster_labels.$BASIS"
  mkdir -p "$r/labels"
  python3 - "$r/results/tables/03_clustering/markers/$f.template.csv" "$r/labels/$f.csv" <<'EOF'
import csv, sys
rows = list(csv.DictReader(open(sys.argv[1])))
with open(sys.argv[2], "w", newline="") as out:
    w = csv.writer(out)
    w.writerow(["cluster_id", "label"])
    for row in rows:
        w.writerow([row["cluster_id"], row["best_guess"] or "Unassigned"])
EOF
}

for id in "${MAKE[@]}"; do
  extra=()
  [ "$id" = smoke_new ] && extra=("${IMAGE[@]}")
  run "$id" discover discover ${extra[@]+"${extra[@]}"} || continue
  draft_labels "$id" || { FAILED+=("$id: labels"); continue; }
  # The cci launch builds the annotate and spatial objects first
  run "$id" cci "annotate, spatial, cci" ${extra[@]+"${extra[@]}"}
  run "$id" de de ${extra[@]+"${extra[@]}"}
done

# Every CSV under results/tables, compared byte for byte
compare_tables() {
  local a="$RUNS/smoke_old/results/tables" b="$RUNS/smoke_new/results/tables" f
  local same=0 differ=() single=()
  while IFS= read -r f; do
    if [ ! -f "$b/$f" ]; then single+=("smoke_old only: $f")
    elif cmp -s "$a/$f" "$b/$f"; then same=$((same + 1))
    else differ+=("$f"); fi
  done < <(cd "$a" 2>/dev/null && find . -name '*.csv' | sed 's|^\./||' | sort)
  while IFS= read -r f; do
    [ -f "$a/$f" ] || single+=("smoke_new only: $f")
  done < <(cd "$b" 2>/dev/null && find . -name '*.csv' | sed 's|^\./||' | sort)
  echo "Tables identical: $same"
  echo "Tables that differ: ${#differ[@]}"
  [ ${#differ[@]} -eq 0 ] || printf '  %s\n' "${differ[@]}"
  echo "Tables in one run only: ${#single[@]}"
  [ ${#single[@]} -eq 0 ] || printf '  %s\n' "${single[@]}"
}

{
  echo "Image: $SIF"
  echo "Runs made now: ${MAKE[*]}"
  echo "Failed launches: ${#FAILED[@]}"
  [ ${#FAILED[@]} -eq 0 ] || printf '  %s\n' "${FAILED[@]}"
  echo "Tasks that ran in the image:"
  for id in smoke_old smoke_new; do
    n=$(find "$SCRATCH/$id" -name .command.run 2>/dev/null | wc -l)
    k=$(find "$SCRATCH/$id" -name .command.run -exec grep -l -F "$SIF" {} + 2>/dev/null | wc -l)
    echo "  $id: $k of $n"
  done
  compare_tables
} | tee "$RUNS/smoke_compare.txt"

#!/usr/bin/env bash
# Checks an image against R/renv.lock, pandoc and Giotto's Python, and loads giotto_dir
#   Environment/check_image.sh <image.sif> [giotto_dir]
set -euo pipefail

[ -f "${1:-}" ] || { echo "usage: Environment/check_image.sh <image.sif> [giotto_dir]" >&2; exit 1; }
SIF="$(cd "$(dirname "$1")" && pwd)/$(basename "$1")"
GDIR=""
if [ -n "${2:-}" ]; then GDIR=$(cd "$2" && pwd); fi
# Apptainer binds the current folder, so the scripts are read from here
cd "$(dirname "$0")/.."

if apptainer exec "$SIF" cat /opt/cosmx-nf/renv.lock | cmp -s - Environment/R/renv.lock; then
  echo "lock    : the image was built from this checkout's Environment/R/renv.lock"
else
  echo "lock    : the image was built from another renv.lock; checking it against this checkout's"
fi
# R with only the image's own library, as in the pipeline's tasks
apptainer exec --cleanenv --env R_LIBS=/none --env R_LIBS_USER=/none "$SIF" \
  Rscript Environment/R/check_image.R Environment/R/renv.lock "$GDIR"

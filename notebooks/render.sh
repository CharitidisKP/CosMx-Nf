#!/usr/bin/env bash
# Knit manual_annotation.Rmd inside the pipeline's container, where Giotto lives.
#
#   ./notebooks/render.sh <run_id> [cluster_basis] [giotto_dir]
# cluster_basis defaults to leiden_clus_res0.3
# giotto_dir defaults to the run's discover object, the subcluster one for sub_*
#
#   ./notebooks/render.sh my_run
#   ./notebooks/render.sh my_run sub_res0.3
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
RUN=${1:?usage: render.sh <run_id> [cluster_basis] [giotto_dir]}
BASIS=${2:-leiden_clus_res0.3}
GDIR=${3:-}
[ -n "$GDIR" ] && GDIR=$(cd "$GDIR" && pwd)   # params are read from inside the container
OUT="manual_annotation.${RUN}.${BASIS}.html"
# Runs live beside the pipeline, or in the project cosmx-nf names in COSMX_NF_RUNS
RUNS=${COSMX_NF_RUNS:-$PIPE/../Runs}

# The image and settings the stages use, with the runs folder passed in
APPTAINERENV_COSMX_NF_RUNS="$RUNS" "$PIPE/Environment/in_image.sh" Rscript -e "
rmarkdown::render('$PIPE/notebooks/manual_annotation.Rmd',
  params = list(run_id = '$RUN', cluster_basis = '$BASIS', giotto_dir = '$GDIR'),
  output_file = '$OUT', output_format = 'html_document')"

echo "wrote $PIPE/notebooks/$OUT"

#!/usr/bin/env bash
# Knit manual_annotation.Rmd inside the pipeline's container, where Giotto lives.
#
#   ./notebooks/render.sh <run_id> [cluster_basis] [giotto_dir]
# defaults: cluster_basis = leiden_clus_res0.3, giotto_dir = the run's discover object
# (objects/subcluster for a sub_* basis)
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

# Same image, libraries and binds as the pipeline, read from conf/site.yaml
site() { sed -n "s/^$1:[[:space:]]*//p" "$PIPE/conf/site.yaml" 2>/dev/null | head -n 1 | tr -d "\"'"; }
SIF=$(site sif)
[ -n "$SIF" ] || { echo "no sif in $PIPE/conf/site.yaml (see conf/site.example.yaml)" >&2; exit 1; }
OPTS=(--cleanenv)
[ -n "$(site r_libs_user)" ]      && OPTS+=(--env "R_LIBS_USER=$(site r_libs_user)")
[ -n "$(site r_libs)" ]           && OPTS+=(--env "R_LIBS=$(site r_libs)")
[ -n "$(site container_python)" ] && OPTS+=(--env "RETICULATE_PYTHON=$(site container_python)")
IFS=, read -r -a BINDS <<<"$(site bind)"
for b in ${BINDS[@]+"${BINDS[@]}"}; do OPTS+=(-B "$b"); done

apptainer exec "${OPTS[@]}" "$SIF" Rscript -e "
rmarkdown::render('$PIPE/notebooks/manual_annotation.Rmd',
  params = list(run_id = '$RUN', cluster_basis = '$BASIS', giotto_dir = '$GDIR'),
  output_file = '$OUT', output_format = 'html_document')"

echo "wrote $PIPE/notebooks/$OUT"

#!/usr/bin/env bash
# Runs a command in the pipeline's image with the settings and folders the stages get
#   Environment/in_image.sh <command> [args...]
#   Environment/in_image.sh Rscript -e 'library(Giotto)'
# Each setting comes from a project's cosmx-nf.yaml, then conf/site.yaml, then conf/server.config
set -euo pipefail
PIPE="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
[ $# -gt 0 ] || { echo "usage: Environment/in_image.sh <command> [args...]" >&2; exit 1; }
U=${USER:-$(id -un)}

# Prints the default of setting $1 in conf/server.config, with the home folder and user filled in
server_default() {
  awk -v k="$1" '$1 == k && $2 == "=" { v = $0; sub(/^[^"]*"/, "", v); sub(/"[^"]*$/, "", v); print v; exit }' \
    "$PIPE/conf/server.config" |
    sed -e "s|\${params.home}|$HOME|g" -e "s|\${System.getenv('USER')}|$U|g"
}

# Prints setting $1 from the first settings file that has it, else its server.config default
setting() {
  local f
  for f in "${COSMX_NF_PROJECT_YAML:-}" "$PIPE/conf/site.yaml"; do
    if [ -n "$f" ] && [ -f "$f" ] && grep -q "^$1:" "$f"; then
      awk -v k="$1" -v q="'" '$0 ~ "^" k ":" { sub("^" k ":[ \t]*", ""); sub(/[ \t]*(#.*)?$/, "")
        gsub(/"/, ""); gsub(q, ""); print; exit }' "$f"
      return 0
    fi
  done
  server_default "$1"
}

SIF=$(setting sif)
PY=$(setting container_python)
LIBS=$(setting r_libs)
LIBS_USER=$(setting r_libs_user)
RTMP=$(setting tmpdir)
BIND=$(setting bind)
[ -f "$SIF" ] || {
  echo "in_image.sh: no image at '$SIF'. cosmx-nf setup pulls it and names it in conf/site.yaml" >&2
  exit 1
}
# Makes the R temp folder the image is given
mkdir -p "$RTMP"

# One -B per comma-separated folder, as in conf/server.config
BINDS=()
FOLDERS=()
IFS=, read -r -a FOLDERS <<<"$BIND"
for b in ${FOLDERS[@]+"${FOLDERS[@]}"}; do BINDS+=(-B "$b"); done

exec apptainer exec --cleanenv \
  --env RETICULATE_PYTHON="$PY" \
  --env R_LIBS="$LIBS" \
  --env R_LIBS_USER="$LIBS_USER" \
  --env TMPDIR="$RTMP" \
  --env VROOM_TEMP_PATH="$RTMP" \
  ${BINDS[@]+"${BINDS[@]}"} \
  "$SIF" "$@"

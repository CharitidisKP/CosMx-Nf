#!/usr/bin/env bash
# Builds the CosMx-Nf image from apptainer/cosmx-nf.def
#   Environment/build.sh [output.sif]        default ./cosmx-nf.sif
set -euo pipefail

out=${1:-cosmx-nf.sif}
case $out in /*) ;; *) out=$PWD/$out ;; esac
[ ! -e "$out" ] || { echo "build.sh: $out exists; give a new path" >&2; exit 1; }

# The recipe's %files paths are relative to this folder
cd "$(dirname "$0")"
[ -s R/renv.lock ] || { echo "build.sh: R/renv.lock is missing, see README.md" >&2; exit 1; }

opts=()
[ "$(id -u)" -eq 0 ] || opts+=(--fakeroot)
apptainer build ${opts[@]+"${opts[@]}"} "$out" apptainer/cosmx-nf.def

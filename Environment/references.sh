#!/usr/bin/env bash
# Downloads or builds the NicheNet networks and the HCA_Kidney and KPMP profiles into assets/
#   Environment/references.sh <cosmx-nf.sif>
set -euo pipefail

sif=${1:?usage: references.sh <cosmx-nf.sif>}
case $sif in /*) ;; *) sif=$PWD/$sif ;; esac
cd "$(dirname "$0")/.."
mkdir -p assets/nichenet assets/profiles
tmp=$(mktemp -d)
trap 'rm -rf "$tmp"' EXIT

# check <file> <md5 | size:bytes>
check() {
  case $2 in
    size:*) [ "$(stat -c %s "$1")" = "${2#size:}" ] ;;
    *)      [ "$(md5sum "$1" | cut -d' ' -f1)" = "$2" ] ;;
  esac
}

# fetch <url> <destination> <md5 | size:bytes>
fetch() {
  if [ -e "$2" ]; then
    if check "$2" "$3"; then echo "have $2"; else
      echo "keep $2: differs from the published file ($3); delete it to download" >&2; fi
    return
  fi
  echo "get  $2"
  curl -fL --retry 3 -o "$2.part" "$1"
  check "$2.part" "$3" || { echo "check failed for $2 ($3)" >&2; rm -f "$2.part"; exit 1; }
  mv "$2.part" "$2"
}

z=https://zenodo.org/records/7074291/files
for f in ligand_target_matrix_nsga2r_final.rds:b09606b04b2d4490418d9028c0e58b9f \
         lr_network_human_21122021.rds:2a155f81e9ffffd5d5e709fe66bcc465 \
         weighted_networks_nsga2r_final.rds:80014dee22df42e98d1c608731c685b5; do
  fetch "$z/${f%%:*}?download=1" "assets/nichenet/${f%%:*}" "${f##*:}"
done

# HCA_Kidney: CellProfileLibrary's .RData at a fixed commit, written out as CSV
if [ -e assets/profiles/HCA_Kidney.csv ]; then
  echo "have assets/profiles/HCA_Kidney.csv"
else
  cpl=https://raw.githubusercontent.com/Nanostring-Biostats/CellProfileLibrary
  cpl=$cpl/9698cc40641940581eac92a2bf21cdcbadb49367/Human/Adult/Kidney_HCA.RData
  fetch "$cpl" "$tmp/Kidney_HCA.RData" 51538d75926ed0cdda23d4070fbc16a4
  apptainer exec "$sif" Rscript -e \
    'e <- new.env(); load(commandArgs(TRUE)[1], envir = e); write.csv(e$profile_matrix, commandArgs(TRUE)[2])' \
    "$tmp/Kidney_HCA.RData" assets/profiles/HCA_Kidney.csv
  echo "made assets/profiles/HCA_Kidney.csv"
fi

# KPMP: the v1.5 atlas from CZ CELLxGENE (1.8 GB), reduced to the profile, then deleted
if [ -e assets/profiles/KPMP.rds ]; then
  echo "have assets/profiles/KPMP.rds"
else
  h5=https://datasets.cellxgene.cziscience.com/08648660-9898-4644-a9a1-ff3a159a0041.h5ad
  fetch "$h5" "$tmp/KPMP_sc_v1.5.h5ad" size:1821593028
  apptainer exec "$sif" Rscript Environment/kpmp_profile.R \
    "$tmp/KPMP_sc_v1.5.h5ad" assets/profiles/KPMP.rds
fi

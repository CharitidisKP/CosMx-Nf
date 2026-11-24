#!/usr/bin/env bash
# Checks Environment/in_image.sh and the helpers that run through it, against a stand-in apptainer
#   tests/test_in_image.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
T=$(cd "$(mktemp -d "${TMPDIR:-/tmp}/in-image-test.XXXXXX")" && pwd -P)
P=$T/Pipeline
export HOME=$T/home USER=cosmxtest$$
trap 'rm -rf "$T" "/tmp/$USER-nf"' EXIT
mkdir -p "$P/conf" "$P/Environment" "$P/notebooks" "$P/shiny" "$T/bin" "$HOME" "$T/images"
cp "$REPO/Environment/in_image.sh" "$P/Environment/"
cp "$REPO/notebooks/render.sh" "$REPO/notebooks/positron_kernel.sh" "$REPO/notebooks/kernel.Rprofile" \
  "$P/notebooks/"
cp "$REPO/shiny/export_app_bundle.sh" "$P/shiny/"
cp "$REPO/conf/server.config" "$P/conf/"
echo image >"$T/images/cosmx-nf.sif"

# The stand-in apptainer writes its arguments one per line, and the variables it would pass in
cat >"$T/bin/apptainer" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" >"$ARGS_OUT"
env | grep '^APPTAINERENV_' | sort >"$ENV_OUT"
EOF
chmod +x "$T/bin/apptainer"
export PATH="$T/bin:$PATH" ARGS_OUT=$T/args ENV_OUT=$T/env
unset COSMX_NF_PROJECT_YAML COSMX_NF_RUNS
cd "$T"

FAILS=0
# check <what> <command...>
check() {
  local what=$1
  shift
  if "$@"; then
    echo "ok    $what"
  else
    echo "FAIL  $what"
    FAILS=$((FAILS + 1))
  fi
}
# run <command...>: runs it with a fresh record of the apptainer call
run() {
  rm -f "$T/args" "$T/env"
  "$@" >"$T/out" 2>&1
  echo $? >"$T/status"
}
# $2 follows $1 somewhere among the arguments apptainer got
pair() {
  awk -v a="$1" -v b="$2" 'prev == a && $0 == b { found = 1 } { prev = $0 } END { exit !found }' \
    "$T/args" 2>/dev/null
}
got() { grep -qxF -- "$1" "$T/args" 2>/dev/null; }
not_got() { ! got "$1"; }
said() { grep -qF -- "$1" "$T/out"; }
failed() { [ "$(cat "$T/status")" != 0 ]; }
passed_in() { grep -qxF -- "$1" "$T/env" 2>/dev/null; }

echo "== settings from conf/site.yaml"
cat >"$P/conf/site.yaml" <<EOF
profile: server
sif: $T/images/cosmx-nf.sif
container_python: /opt/py/bin/python
python_path: /opt/py/bin/python
r_libs: /none
r_libs_user: /none
tmpdir: $T/rtmp
# bind: /data,/scratch
EOF
run "$P/Environment/in_image.sh" Rscript -e 1
check "it runs apptainer exec" test "$(head -n 1 "$T/args")" = exec
check "with --cleanenv" got --cleanenv
check "on the image site.yaml names" pair "$T/images/cosmx-nf.sif" Rscript
check "the command and its arguments come last" test "$(tail -n 3 "$T/args" | tr '\n' ' ')" = "Rscript -e 1 "
check "with the image's Python" pair --env RETICULATE_PYTHON=/opt/py/bin/python
check "no host R library" pair --env R_LIBS=/none
check "in either R library variable" pair --env R_LIBS_USER=/none
check "the R temp folder from the settings" pair --env "TMPDIR=$T/rtmp"
check "for vroom too" pair --env "VROOM_TEMP_PATH=$T/rtmp"
check "which it makes" test -d "$T/rtmp"
check "a commented bind binds nothing" not_got -B

echo "== bind"
printf 'sif: %s\nbind: /data,/scratch\n' "$T/images/cosmx-nf.sif" >"$P/conf/site.yaml"
run "$P/Environment/in_image.sh" true
check "each folder in bind is bound" pair -B /data
check "comma-separated" pair -B /scratch

echo "== defaults from conf/server.config"
check "a key no settings file has is passed empty" pair --env RETICULATE_PYTHON=
check "the R temp folder is the server.config default, with the user filled in" \
  pair --env "TMPDIR=/tmp/$USER-nf/Rtmp"
printf 'sif: %s\nbind:\n' "$T/images/cosmx-nf.sif" >"$P/conf/site.yaml"
run "$P/Environment/in_image.sh" true
check "an empty bind binds nothing" not_got -B
rm "$P/conf/site.yaml"
run "$P/Environment/in_image.sh" true
check "without site.yaml there is no image" said "no image at ''"
check "a missing image stops it" failed
check "before apptainer runs" test ! -e "$T/args"
run "$P/Environment/in_image.sh"
check "no command is refused" said "usage: Environment/in_image.sh"

echo "== quotes, comments and a project's settings"
cat >"$P/conf/site.yaml" <<EOF
sif: "$T/images/cosmx-nf.sif"   # the image
container_python: '/opt/py/bin/python'
EOF
run "$P/Environment/in_image.sh" true
check "quotes and a trailing comment are dropped" pair "$T/images/cosmx-nf.sif" true
check "single quotes too" pair --env RETICULATE_PYTHON=/opt/py/bin/python
echo image >"$T/images/project.sif"
printf 'sif: %s\ninput: samplesheet.csv\n' "$T/images/project.sif" >"$T/cosmx-nf.yaml"
COSMX_NF_PROJECT_YAML=$T/cosmx-nf.yaml run "$P/Environment/in_image.sh" true
check "a project's setting wins" pair "$T/images/project.sif" true
check "and the rest come from site.yaml" pair --env RETICULATE_PYTHON=/opt/py/bin/python

echo "== the helpers"
cat >"$P/conf/site.yaml" <<EOF
sif: $T/images/cosmx-nf.sif
container_python: /opt/py/bin/python
r_libs: /none
r_libs_user: /none
tmpdir: $T/rtmp
EOF
run "$P/notebooks/render.sh" r1 sub_res0.3
check "render.sh knits in the image" pair "$T/images/cosmx-nf.sif" Rscript
check "the annotation notebook" grep -qF "rmarkdown::render('$P/notebooks/manual_annotation.Rmd'," "$T/args"
check "for the run and clustering given" grep -qF "params = list(run_id = 'r1', cluster_basis = 'sub_res0.3'" "$T/args"
check "with the runs folder passed in" passed_in "APPTAINERENV_COSMX_NF_RUNS=$P/../Runs"
mkdir -p "$T/Runs/r1/results/objects/discover/giotto"
COSMX_NF_RUNS=$T/Runs run "$P/shiny/export_app_bundle.sh" r1
check "export_app_bundle.sh exports in the image" pair "$T/images/cosmx-nf.sif" Rscript
check "with export_bundle.R" pair --vanilla "$P/bin/export_bundle.R"
check "from the run's object" pair --input "$T/Runs/r1/results/objects/discover/giotto"
check "with the image's Python" pair --env RETICULATE_PYTHON=/opt/py/bin/python
mkdir -p "$T/positron/resources/scripts" "$T/positron/resources/ark"
touch "$T/positron/resources/scripts/startup.R"
printf '#!/bin/sh\n' >"$T/positron/resources/ark/ark"
chmod +x "$T/positron/resources/ark/ark"
POSITRON=1 run "$P/notebooks/positron_kernel.sh" --connection_file c.json \
  --startup-file "$T/positron/resources/scripts/startup.R"
check "the kernel runs Positron's ark in the image" pair "$T/images/cosmx-nf.sif" sh
check "after R's library paths" pair -c '. "$R_HOME/etc/ldpaths" && exec "$@"'
check "ark itself" pair ark "$T/positron/resources/ark/ark"
check "with Positron's arguments" pair --startup-file "$T/positron/resources/scripts/startup.R"
check "Positron's variables pass through" passed_in APPTAINERENV_POSITRON=1
check "and the kernel's R profile" passed_in "APPTAINERENV_R_PROFILE_USER=$P/notebooks/kernel.Rprofile"
run "$P/notebooks/positron_kernel.sh" --startup-file "$T/nowhere/scripts/startup.R"
check "no ark beside --startup-file is refused" said "no ark beside --startup-file"
check "before apptainer runs" test ! -e "$T/args"

echo
if [ "$FAILS" -eq 0 ]; then
  echo "all passed"
else
  echo "$FAILS failed"
  exit 1
fi

#!/usr/bin/env bash
# Checks cosmx-nf against stand-ins for nextflow, apptainer and the helper scripts
#   tests/test_cosmx-nf.sh
set -uo pipefail

REPO=$(cd "$(dirname "$0")/.." && pwd)
T=$(cd "$(mktemp -d "${TMPDIR:-/tmp}/cosmx-nf-test.XXXXXX")" && pwd -P)
trap 'rm -rf "$T"' EXIT
P=$T/Pipeline
R=$T/Runs
mkdir -p "$P/conf" "$P/notebooks" "$P/shiny" "$P/Environment" "$P/assets" "$P/.nextflow" \
  "$T/bin" "$T/home/bin"
cp "$REPO/cosmx-nf" "$REPO/nextflow.config" "$P/"
cp "$REPO"/conf/*.config "$P/conf/"
cp "$REPO/assets/no_label_overrides.csv" "$P/assets/"
printf 'sample_id,include,data_dir\nS1,TRUE,/data/S1\n' >"$P/assets/samplesheet.csv"

# The stand-in nextflow writes its arguments one per line
cat >"$T/bin/nextflow" <<'EOF'
#!/usr/bin/env bash
printf '%s\n' "$@" >"$ARGS_OUT"
EOF
# The stand-in render.sh writes a draft where the notebook does
cat >"$P/notebooks/render.sh" <<'EOF'
#!/usr/bin/env bash
d="$COSMX_NF_RUNS/$1/results/tables/06_annotation/labels"
mkdir -p "$d"
printf 'cluster_id,label\n1,PT\n' >"$d/cluster_labels.$2.csv"
EOF
# The stand-in bundle export and references.sh note how they were called
cat >"$P/shiny/export_app_bundle.sh" <<'EOF'
#!/usr/bin/env bash
echo "$COSMX_NF_RUNS|$*" >"$CALLED"
echo "project settings: ${COSMX_NF_PROJECT_YAML:-none}" >>"$CALLED"
EOF
cat >"$P/Environment/references.sh" <<'EOF'
#!/usr/bin/env bash
echo "references $*" >>"$CALLED"
EOF
# The stand-in apptainer pulls by writing the file
cat >"$T/bin/apptainer" <<'EOF'
#!/usr/bin/env bash
echo "apptainer $*" >>"$CALLED"
if [ "$1" = pull ]; then echo image >"$2"; fi
EOF
chmod +x "$T/bin/nextflow" "$T/bin/apptainer" "$P/notebooks/render.sh" \
  "$P/shiny/export_app_bundle.sh" "$P/Environment/references.sh"
export NXF_BIN=$T/bin/nextflow NXF_SCRATCH_ROOT=$T/scratch ARGS_OUT=$T/args CALLED=$T/called
# Outside any project folder the tests run from
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
launch() {
  rm -f "$T/args" "$T/called"
  "$P/cosmx-nf" "$@" >"$T/out" 2>&1
}
# $2 follows $1 among the arguments nextflow got
passed() { grep -A1 -xF -- "$1" "$T/args" 2>/dev/null | tail -n 1 | grep -qxF -- "$2"; }
not_passed() { ! grep -qxF -- "$1" "$T/args" 2>/dev/null; }
said() { grep -qF -- "$1" "$T/out"; }
not_said() { ! said "$1"; }
called() { grep -qF -- "$1" "$T/called" 2>/dev/null; }
not_called() { ! called "$1"; }
# history <file> <session> <status> <command tail>: a line of Nextflow's history
history() {
  printf '2026-10-30 10:00:00\t1m 0s\tname_%s\t%s\t-\t%s\tnextflow run main.nf %s\n' \
    "${2:0:4}" "$3" "$2" "$4" >>"$1"
}
H=$P/.nextflow/history
U1=11111111-1111-1111-1111-111111111111
U2=22222222-2222-2222-2222-222222222222
U3=33333333-3333-3333-3333-333333333333
U4=44444444-4444-4444-4444-444444444444

echo "== stage word, profile, first launch"
launch r1 discover
check "the stage word becomes --stage" passed --stage discover
check "the profile defaults to server" passed -profile server
check "the run_id is set" passed --run_id r1
check "results go beside the pipeline" passed --outdir "$P/../Runs/r1/results"
check "a first launch starts fresh" not_passed -resume
check "it says why" said "no earlier discover launch for r1"
launch r1 --stage spatial
check "--stage still works" passed --stage spatial

echo "== resuming"
# A discover launch from before the stage was always passed carries no --stage
history "$H" "$U1" OK "--outdir x --run_id r1 -profile server"
history "$H" "$U3" OK "--outdir x --run_id r10 -profile server --stage discover"
launch r1 discover
check "a stage resumes its last launch for this run" passed -resume "$U1"
launch r1 discover -resume
check "a bare -resume is the same" passed -resume "$U1"
launch r1 discover --fresh
check "--fresh starts over" not_passed -resume
launch r1 discover --fresh -resume
check "--fresh with -resume is refused" said "contradict each other"
launch r1 discover -resume "$U2"
check "an explicit session is kept" passed -resume "$U2"
mkdir -p "$R/r1/results/tables/08_cci" "$R/r1/results/objects/spatial/giotto"
touch -t 202610300900 "$R/r1/results/tables/08_cci/liana.csv"
touch -t 202610301000 "$R/r1/results/objects/spatial/giotto/gobject.RDS"
history "$H" "$U2" OK "--outdir x --run_id r1 -profile server --stage cci"
launch r1 cci
check "an input object newer than the outputs starts fresh" not_passed -resume
check "it says why" said "objects/spatial/giotto is newer than every cci output"
touch -t 202610301100 "$R/r1/results/tables/08_cci/liana.csv"
launch r1 cci
check "otherwise cci resumes its own last launch" passed -resume "$U2"

echo "== parameters"
launch r1 discover --input sheet.csv --seed 7
check "declared parameters pass" passed --input sheet.csv
check "every one of them" passed --seed 7
launch r1 discover --no_such_param 1
check "an undeclared parameter is refused" said "unknown parameter: --no_such_param"
launch r1 spatail
check "a mistyped stage is refused" said "no stage or command 'spatail'"
launch r1 -profile server
check "options without a stage are refused" said "name a stage"
launch "a/b" discover
check "a run_id with a slash is refused" said "a run_id is"

echo "== profile and site.yaml"
printf 'profile: "slurm"   # the cluster\nprofile_dir: /x\n' >"$P/conf/site.yaml"
launch r1 discover
check "the profile comes from the settings" passed -profile slurm
check "conf/site.yaml is the params file" passed -params-file "$P/conf/site.yaml"
launch r1 discover -profile server
check "-profile overrides the settings" passed -profile server
rm "$P/conf/site.yaml"

echo "== nextflow off PATH"
cp "$T/bin/nextflow" "$T/home/bin/nextflow"
(
  unset NXF_BIN
  export HOME=$T/home PATH=/usr/bin:/bin
  launch r1 discover
)
check "nextflow is found in ~/bin" passed --stage discover
check "it says where" said "$T/home/bin/nextflow (not on PATH)"

echo "== status"
launch r9
check "a new run's next command is discover" said "next    : cosmx-nf r9 discover"
O=$R/r2/results/objects
L=$R/r2/labels
mkdir -p "$O/discover/giotto" "$O/annotate/giotto" "$O/spatial/giotto" "$L" "$R/r2/results/tables/08_cci"
touch -t 202610010900 "$O/discover/giotto/gobject.RDS"
touch -t 202610020900 "$O/annotate/giotto/gobject.RDS"
touch -t 202610030900 "$O/spatial/giotto/gobject.RDS"
launch r2
check "without labels the next command drafts them" said "next    : cosmx-nf r2 labels <clustering>"
mkdir -p "$R/r2/results/tables/03_clustering/checks"
echo leiden_clus_res0.4 >"$R/r2/results/tables/03_clustering/checks/picked_resolution.txt"
launch r2
check "with a pick the next command drafts it" said "next    : cosmx-nf r2 labels (drafts leiden_clus_res0.4"
check "without a site status names none" not_said "report  :"
mkdir -p "$R/r2/results/report"
touch "$R/r2/results/report/index.html"
launch r2
check "status names the HTML report" said "report  : $R/r2/results/report/index.html"
printf 'cluster_id,label\n1,PT\n2,TAL\n' >"$L/cluster_labels.leiden_clus_res0.3.csv"
printf 'cell_ID,label,confidence,note\nc1,Podocyte,low,checked\n' >"$L/label_overrides.csv"
# main.nf's digest of these labels: md5 of "<clustering>\n<labels>\n<edits>"
echo 7f2e6431922f3a99a226014ac53b08bb >"$O/annotate/label_state.txt"
launch r2
check "status finds the labels annotate was built from" said "annotate is up to date with them"
check "outputs without a launch record count as done" said "(launched from another folder)"
check "the next command is the stage not yet run" said "next    : cosmx-nf r2 de"
mkdir -p "$R/r2/results/tables/09_de"
touch -t 202610020900 "$R/r2/results/tables/09_de/smide.csv"
launch r2
check "results older than the object they read are outdated" said "outdated: finished before"
check "and that stage runs next" said "next    : cosmx-nf r2 de"
printf '3,DCT\n' >>"$L/cluster_labels.leiden_clus_res0.3.csv"
launch r2
check "status sees labels edited after annotate" said "edited since annotate ran"
check "and sends them through cci" said "next    : cosmx-nf r2 cci (annotate and spatial run first)"
printf 'cluster_id,label\n1,PT\n' >"$L/cluster_labels.sub_res0.3.csv"
launch r2
check "two labelled clusterings need cluster_basis.txt" said "write the one cell_type comes from"

echo "== labels and bundle"
mkdir -p "$R/r3/results/objects/discover/giotto"
touch "$R/r3/results/objects/discover/giotto/gobject.RDS"
launch r3 labels
check "labels copies the draft into labels/" test -f "$R/r3/labels/cluster_labels.leiden_clus_res0.3.csv"
printf 'cluster_id,label\n1,Curated\n' >"$R/r3/labels/cluster_labels.leiden_clus_res0.3.csv"
launch r3 labels
check "labels keeps curated labels" grep -q Curated "$R/r3/labels/cluster_labels.leiden_clus_res0.3.csv"
check "and says where the new draft is" said "kept the curated"
launch r3 labels sub_res0.3
check "a subclustering needs the subcluster object" said "r3 has no subcluster object yet"
mkdir -p "$R/r4/results/objects/discover/giotto" "$R/r4/results/tables/03_clustering/checks"
touch "$R/r4/results/objects/discover/giotto/gobject.RDS"
echo leiden_clus_res0.4 >"$R/r4/results/tables/03_clustering/checks/picked_resolution.txt"
launch r4 labels
check "labels without a clustering drafts the pick" test -f "$R/r4/labels/cluster_labels.leiden_clus_res0.4.csv"
check "and says it is the pick" said "the stability check's pick"
launch r3 bundle all xz objects/annotate
check "bundle runs the export with the run's folder" called "$P/../Runs|r3 all xz objects/annotate"
check "outside a project the helpers get no project settings" called "project settings: none"
COSMX_NF_PROJECT_YAML=/stale/cosmx-nf.yaml launch r3 bundle
check "even one left in the shell" called "project settings: none"

echo "== install"
(
  export HOME=$T/home PATH=/usr/bin:/bin
  launch install
)
check "install links cosmx-nf into ~/.local/bin" test -L "$T/home/.local/bin/cosmx-nf"
check "and says how to put it on PATH" said "export PATH=\"$T/home/.local/bin:\$PATH\""
rm -f "$T/args"
"$T/home/.local/bin/cosmx-nf" r1 discover >"$T/out" 2>&1
check "the link runs the pipeline it points to" passed --outdir "$P/../Runs/r1/results"

echo "== init and projects"
launch init "$T/study"
check "init writes the project's settings" grep -qx "input: samplesheet.csv" "$T/study/cosmx-nf.yaml"
check "a samplesheet with the columns only" test "$(cat "$T/study/samplesheet.csv")" = "sample_id,include,data_dir"
check "and a Runs folder" test -d "$T/study/Runs"
launch init "$T/study"
check "init refuses a project twice" said "is a project already"
mkdir -p "$T/study/sub"
rm -f "$T/args"
(cd "$T/study/sub" && "$P/cosmx-nf" p1 discover >"$T/out" 2>&1)
check "in a project, results go to its Runs" passed --outdir "$T/study/Runs/p1/results"
check "the project's settings are the params file" passed -params-file "$T/study/cosmx-nf.yaml"
check "Nextflow runs the pipeline's main.nf" passed run "$P/main.nf"
printf 'profile: server\nsif: /img/cosmx-nf.sif\ninput: /x/real.csv\n' >"$P/conf/site.yaml"
rm -f "$T/args"
(cd "$T/study" && "$P/cosmx-nf" p1 discover >"$T/out" 2>&1)
S=$T/study/Runs/p1/logs/settings.yaml
check "with a site.yaml, the two are merged" passed -params-file "$S"
check "the machine's keys are kept" grep -qx "sif: /img/cosmx-nf.sif" "$S"
check "the project's keys win" grep -qx "input: samplesheet.csv" "$S"
check "each key once" test "$(grep -c '^input:' "$S")" = 1
history "$H" "$U4" OK "--outdir x --run_id p1 -profile server --stage discover"
rm -f "$T/args"
(cd "$T/study" && "$P/cosmx-nf" p1 discover >"$T/out" 2>&1)
check "a project reads its own history, not the pipeline folder's" not_passed -resume
mkdir -p "$T/study/.nextflow"
history "$T/study/.nextflow/history" "$U4" OK "--outdir x --run_id p1 -profile server --stage discover"
rm -f "$T/args"
(cd "$T/study" && "$P/cosmx-nf" p1 discover >"$T/out" 2>&1)
check "and resumes from it" passed -resume "$U4"
rm -f "$T/called"
(cd "$T/study" && "$P/cosmx-nf" p1 bundle >"$T/out" 2>&1)
check "bundle gets the project's runs folder" called "$T/study/Runs|p1"
check "and the project's settings for the image" called "project settings: $T/study/cosmx-nf.yaml"
rm "$P/conf/site.yaml"

echo "== setup"
(
  export PATH=/usr/bin:/bin
  launch setup
)
check "setup needs apptainer" said "apptainer is not installed here"
(
  export PATH=$T/bin:/usr/bin:/bin
  launch setup
)
check "setup pulls the image for this version" called "apptainer pull $P/images/cosmx-nf_0.1.0.sif oras://ghcr.io/charitidiskp/cosmx-nf:0.1.0"
check "fetches the references with it" called "references $P/images/cosmx-nf_0.1.0.sif"
check "and writes conf/site.yaml" grep -qx "sif: $P/images/cosmx-nf_0.1.0.sif" "$P/conf/site.yaml"
check "with the image's Python" grep -qx "python_path: /opt/py/bin/python" "$P/conf/site.yaml"
check "and no host R library" grep -qx "r_libs: /none" "$P/conf/site.yaml"
check "for either R variable" grep -qx "r_libs_user: /none" "$P/conf/site.yaml"
echo "profile: slurm" >>"$P/conf/site.yaml"
(
  export PATH=$T/bin:/usr/bin:/bin
  launch setup --no-references
)
check "an image already pulled is kept" said "(already here)"
check "--no-references skips them" not_called references
check "an existing conf/site.yaml is left alone" grep -qx "profile: slurm" "$P/conf/site.yaml"
check "and the lines it needs are shown" said "r_libs_user: /none"
echo image >"$T/own.sif"
(
  export PATH=$T/bin:/usr/bin:/bin
  launch setup --sif "$T/own.sif" --no-references
)
check "--sif uses an image already on disk" said "sif: $T/own.sif"
check "and pulls nothing" test ! -e "$T/called"

echo
if [ "$FAILS" -eq 0 ]; then
  echo "all passed"
else
  echo "$FAILS failed"
  exit 1
fi

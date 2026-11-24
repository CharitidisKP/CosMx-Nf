# Environment

One Apptainer image holds everything the pipeline runs on:

- R 4.5.1.
- Every R package, at the version and source recorded in `R/renv.lock`.
- The Python that Giotto calls through reticulate, at `/opt/py`.

Nothing is installed on the host apart from Apptainer and Nextflow.

## Files

| File | What it is |
|---|---|
| `apptainer/cosmx-nf.def` | Image recipe, built on `rocker/rstudio:4.5.1` pinned by digest |
| `R/renv.lock` | Every R package the pipeline and its notebooks load, with its version and its source: CRAN, a Bioconductor release, or a GitHub commit |
| `R/install.R` | Installs the lock into the image during the build |
| `R/snapshot.R` | Writes `renv.lock` from an environment the pipeline runs in |
| `python/conda-explicit.txt` | The Python env's conda packages, exact builds |
| `python/requirements.txt` | Its pip packages, exact versions |
| `build.sh` | Builds the image |
| `check_image.sh`, `R/check_image.R` | Checks an image: every locked package installed and loading, pandoc, Giotto's Python |
| `references.sh` | Downloads or builds the reference files in `assets/` |
| `kpmp_profile.R` | Builds the KPMP reference profile from the KPMP atlas |

## Get the image

The image is published on GitHub's container registry, one tag per pipeline version (`version`
in `nextflow.config`). `cosmx-nf setup` pulls it, fetches the reference files and writes
`conf/site.yaml`. By hand, without root:

```bash
apptainer pull cosmx-nf_0.1.0.sif oras://ghcr.io/charitidiskp/cosmx-nf:0.1.0
```

## Build the image

Building it yourself needs Linux with Apptainer 1.1 or later, and root or `--fakeroot`:

```bash
Environment/build.sh /path/to/cosmx-nf.sif
apptainer test /path/to/cosmx-nf.sif
```

The test loads the pipeline's main R packages and Giotto's Python modules. Most CRAN packages arrive as binaries. The Bioconductor and GitHub packages compile, so the build takes an hour or more and runs at most 8 installs at a time.

`build.sh` never overwrites an existing file. Apptainer unpacks the build in `$APPTAINER_TMPDIR` (default `/tmp`). That needs about 15 GB free on a local disk. A network home is far slower.

## Check an image

Before using or publishing an image, check it as the pipeline's tasks see it, with no R library from the host:

```bash
Environment/check_image.sh /path/to/cosmx-nf.sif [/path/to/a/run/results/objects/discover/giotto]
```

It reports a problem when R is not the lock's version, or when any package in `R/renv.lock` is missing, at another version, or fails to load. It also reports a problem when rmarkdown cannot render a page (pandoc), or when Giotto's Python fails to import or Giotto rejects it. Given a saveGiotto folder, it also loads that object, as `cosmx-nf <run> labels` and `bundle` do. Loading every package takes a few minutes.

## Test the image against the current setup

```bash
Environment/smoke_test.sh /path/to/cosmx-nf.sif /path/to/assets
```

It runs every stage on `assets/samplesheet.smoke.csv` twice: as `smoke_old` with the server profile as it is, and as `smoke_new` with everything from the image. The smoke sheet holds three FOVs from each of two subjects' baseline and follow-up samples. The labels come from the diagnostics' best guess per cluster. Before any run, it prints the R libraries the image's tasks will see, and stops if one is under your home folder. The second argument is a folder holding the reference `profiles/` and `nichenet/` (default `assets/`).

The runs land in `../Runs/smoke_old` and `../Runs/smoke_new`. `../Runs/smoke_compare.txt` lists any failed launch, how many tasks ran in the image, and the tables that differ between the two runs. It takes a few hours, so run it in tmux.

To redo one run, delete it from `../Runs` and from the scratch folder, then name it after the assets folder:

```bash
Environment/smoke_test.sh /path/to/cosmx-nf.sif /path/to/assets smoke_new
```

The other run is kept as it is and compared again.

## Fetch the reference files

```bash
Environment/references.sh /path/to/cosmx-nf.sif
```

- `assets/nichenet/`: the three human NicheNet v2 networks from [Zenodo record 7074291](https://doi.org/10.5281/zenodo.7074291), used by the `cci` stage.
- `assets/profiles/HCA_Kidney.csv`: `Human/Adult/Kidney_HCA` from [CellProfileLibrary](https://github.com/Nanostring-Biostats/CellProfileLibrary) at a fixed commit, used for supervised cell typing.
- `assets/profiles/KPMP.rds`: built by `kpmp_profile.R` from the KPMP single-cell atlas v1.5 on [CZ CELLxGENE](https://cellxgene.cziscience.com/collections/0f528c8a-a25c-4840-8fa3-d156fa11086f). The atlas is one fixed dataset version (1.8 GB) and is deleted after the build. The profile holds the mean expression per `author_cell_type` for every cell type with at least 200 cells. It is also used for supervised cell typing.

Each download is checked against its published md5. The KPMP atlas has no published md5, so it is checked by size. A file that is already in place is never overwritten. If it fails the check, the script reports it and keeps it.

**NicheNet files downloaded before this script.** The copies the pipeline has used so far differ from the files Zenodo publishes:

- `ligand_target_matrix`: 258,545,184 bytes, against 262,116,605 published.
- `weighted_networks`: 73,689,658 bytes, against 73,662,019 published.
- `lr_network`: matches.

Zenodo has published only one version of the record (13 September 2022). So those copies changed after they were downloaded, possibly by being saved again from R. `references.sh` leaves them in place, but a fresh setup gets the published files. Until the two sets are shown to hold the same networks, `cci` results from the old copies may not reproduce exactly from the published files.

## Publish the image

After a rebuild, raise `version` in `nextflow.config`, then push the image under that tag. The
login asks for a GitHub token with the `write:packages` scope:

```bash
apptainer registry login --username <github-user> oras://ghcr.io
apptainer push cosmx-nf.sif oras://ghcr.io/charitidiskp/cosmx-nf:0.1.0
```

A new package on GitHub starts private. In its Package settings, change the visibility to
public and connect it to the repository, so the pull above needs no login.

## Point the pipeline at the image

`cosmx-nf setup --sif /path/to/cosmx-nf.sif` writes these into `conf/site.yaml`. By hand, set
them there, or in the profile you run with:

```yaml
sif: /path/to/cosmx-nf.sif
python_path: /opt/py/bin/python
container_python: /opt/py/bin/python
r_libs: /none
r_libs_user: /none
```

Every package is inside the image. `/none` keeps out the host R libraries that a profile names,
such as the ones in `conf/server.config`.

## Add or update a package

1. Install it into a library on top of the image:

   ```bash
   mkdir -p ~/Rlibs_extra
   apptainer exec --env R_LIBS=$HOME/Rlibs_extra cosmx-nf.sif \
     Rscript -e 'install.packages("pkg", lib = "~/Rlibs_extra")'
   ```

2. Record the environment with that library in place:

   ```bash
   apptainer exec --cleanenv --env R_LIBS=$HOME/Rlibs_extra cosmx-nf.sif \
     Rscript - < Environment/R/snapshot.R > Environment/R/renv.lock
   ```

   The script reports, on stderr, any package recorded at a different version from the one R loads. It also reports any package whose source renv cannot install again, such as a local tarball.

3. Rebuild the image.

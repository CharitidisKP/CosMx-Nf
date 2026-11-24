// Rebuilds the HTML report from this stage's tables and the results tree
process REPORT {
    tag   { stage }
    label 'process_light'
    // Reruns every launch, and a report error never fails the stage
    cache false
    errorStrategy 'ignore'
    // Prints the build's warnings in the console
    debug true
    publishDir path: { "${params.outdir}" }, mode: 'copy', pattern: "report", overwrite: true

    input:
        val(stage)
        val(commit)
        // One numbered folder per file, so two tables with the same name cannot collide
        path(current, stageAs: 'current*/*')
        path(samplesheet)
        tuple path(site_src, stageAs: 'site_src'), path(config, stageAs: 'pipeline.config'),
              path(no_overrides, stageAs: 'no_label_overrides.csv'), path(fig_style, stageAs: 'cosmx_figures.R'),
              path(label_colours, stageAs: 'label_colours.csv')

    output:
        path "report", optional: true

    script:
    def labels = params.labels_dir ?: "${params.outdir}/../labels"
    def updated = new Date().format('yyyy-MM-dd HH:mm')
    """
    # The site builds from a copy of report/, named after the run
    cp -RL site_src site
    sed -i 's|^  title: .*|  title: "CosMx-Nf run ${params.run_id}"|' site/_quarto.yml
    export COSMX_RESULTS='${params.outdir}' COSMX_CURRENT="\$PWD/current" COSMX_STAGE='${stage}' \\
        COSMX_RUN_ID='${params.run_id}' COSMX_COMMIT='${commit}' COSMX_UPDATED='${updated}' \\
        COSMX_CHOSEN='${params.chosen_resolutions}' COSMX_UNTYPED='${params.untyped_celltypes}' \\
        COSMX_SAMPLESHEET="\$PWD/${samplesheet}" COSMX_CONFIG="\$PWD/pipeline.config" \\
        COSMX_NO_OVERRIDES="\$PWD/no_label_overrides.csv" COSMX_LABELS='${labels}' \\
        COSMX_FIG_STYLE="\$PWD/cosmx_figures.R" COSMX_LABEL_COLOURS="\$PWD/label_colours.csv" \\
        COSMX_HISTORY='${workflow.launchDir}/.nextflow/history' NO_COLOR=1
    # The whole site, then any page it did not reach on its own, the full log kept
    (cd site && quarto render) > site.log 2>&1 || {
        tail -n 40 site.log >&2
        for q in site/*.qmd; do
            p=\$(basename "\$q" .qmd)
            [ -f "site/_site/\$p.html" ] || (cd site && quarto render "\$p.qmd") >> site.log 2>&1 ||
                echo "REPORT: page \$p did not build"
        done
    }
    # Each distinct Quarto error, with the line under a FATAL that names its cause
    sed 's/\\x1b\\[[0-9;]*m//g' site.log | awk '/^FATAL/ { getline n; print n; next } /^ERROR/' |
        awk '!seen[\$0]++' | sed 's|^|REPORT: quarto: |' || true
    if [ -f site/_site/index.html ]; then
        mv site/_site report
        grep -l 'cell-output-error' report/*.html 2>/dev/null | sed 's|^|REPORT: an error shows on |' || true
        # Each section a page could not show, with its error
        awk '/^Not shown: / { getline m; print "REPORT: a section failed on " FILENAME ": " m }' report/*.html || true
    else
        echo "REPORT: the HTML report did not build, results/report keeps the last one"
    fi
    """
}

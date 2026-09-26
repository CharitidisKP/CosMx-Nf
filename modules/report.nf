// Rewrites run_report.md at the end of a stage. This stage's tables arrive through channels,
// because publishing is asynchronous; earlier stages' tables are read from the results tree.
// The update time comes from Nextflow, because the container clock runs on UTC.
process REPORT {
    tag   { stage }
    label 'process_light'
    // Every launch rebuilds the report from the current results tree. A report error must not
    // fail a stage whose results are already published; --stage report rebuilds it
    cache false
    errorStrategy 'ignore'
    publishDir path: { "${params.outdir}" }, mode: 'copy', pattern: "run_report.md", overwrite: true

    input:
        val(stage)
        val(commit)
        // One numbered folder per file, so two tables with the same name cannot collide
        path(current, stageAs: 'current*/*')
        path(samplesheet)

    output:
        path "run_report.md"

    script:
    """
    run_report.R --results '${params.outdir}' --current current --samplesheet '${samplesheet}' \\
        --run_id '${params.run_id}' --stage '${stage}' --commit '${commit}' \\
        --basis '${params.cluster_basis ?: ""}' --untyped '${params.untyped_celltypes}' \\
        --chosen_resolutions '${params.chosen_resolutions}' \\
        --updated '${new Date().format('yyyy-MM-dd HH:mm')}' --outdir .
    """
}

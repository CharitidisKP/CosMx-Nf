process SMIDE_META {
    tag   'smide_meta'
    label 'process_light'
    publishDir path: { "${params.outdir}/figures/09_de" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/09_de" }, mode: 'copy', pattern: "smide_meta*.csv"

    input:
        path de_files

    output:
        path "smide_meta.csv",         emit: meta
        path "smide_meta_summary.csv", emit: summary
        path "*.png", optional: true

    script:
    """
    smide_meta.R --input . --fdr ${params.smide_fdr} --outdir .
    """
}

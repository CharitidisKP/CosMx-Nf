process COMPOSITION {
    tag   'composition'
    label 'process_light'
    publishDir path: { "${params.outdir}/figures/06_annotation" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/06_annotation" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path samplesheet

    output:
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    composition.R --input '${obj}' --samplesheet '${samplesheet}' \\
        --compare_by '${params.compare_by}' --compare_baseline '${params.compare_baseline}' \\
        --compare_within '${params.compare_within}' --pair_by '${params.pair_by}' \\
        --exclude '${params.untyped_celltypes}' \\
        --outdir . --python '${params.python_path}'
    """
}

process QC {
    tag   { meta.id }
    label 'process_qc'
    publishDir path: { "${params.outdir}/tables/01_qc" }, mode: 'copy', pattern: "*.csv"
    publishDir path: { "${params.outdir}/figures/01_qc" }, mode: 'copy', pattern: "*.png"

    input:
        tuple val(meta), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(meta), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    qc.R --sample_id '${meta.id}' --input '${obj}' \\
        --gene_min_cells ${params.gene_min_cells} --cell_min_genes ${params.cell_min_genes} \\
        --count_cap ${params.count_cap} --count_quantile ${params.count_quantile} \\
        --area_max ${params.area_max} --fov_min_count ${params.fov_min_count} \\
        --split_ratio_threshold ${params.split_ratio_threshold} \\
        --outdir . --python '${params.python_path}'
    """
}

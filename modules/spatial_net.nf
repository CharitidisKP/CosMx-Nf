process SPATIAL_NET {
    tag   'spatial_net'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/figures/07_spatial/network" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/07_spatial/network" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    spatial_net.R --input '${obj}' --celltype_column '${params.celltype_column}' \\
        --max_delaunay_dist ${params.max_delaunay_dist} \\
        --prox_sim ${params.prox_sim} --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

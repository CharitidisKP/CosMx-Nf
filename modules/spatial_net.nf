process SPATIAL_NET {
    tag   'spatial_net'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/figures/07_spatial/network" }, mode: 'copy', pattern: "*.{png,pdf,jpg}"
    publishDir path: { "${params.outdir}/tables/07_spatial/network" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "*.{png,pdf,jpg}"

    script:
    """
    spatial_net.R --input '${obj}' --celltype_column '${params.celltype_column}' \\
        --max_delaunay_dist ${params.max_delaunay_dist} \\
        --prox_sim ${params.prox_sim} --seed ${params.seed} \\
        --um_per_px ${params.um_per_px} --immune_types '${params.immune_types}' \\
        --untyped '${params.untyped_celltypes}' --label_colours '${params.label_colours}' \\
        --fig_format '${params.fig_format}' --outdir . --python '${params.python_path}'
    """
}

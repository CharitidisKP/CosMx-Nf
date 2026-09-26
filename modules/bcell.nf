process BCELL {
    tag   'bcell'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/figures/07_spatial/bcell" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/07_spatial/bcell" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        path "*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    bcell.R --input '${obj}' --celltype_column '${params.celltype_column}' --bcell_regex '${params.bcell_regex}' \\
        --min_cells ${params.bcell_min_cells} --markers '${params.markers}' --n_hvgs ${params.n_hvgs} \\
        --n_pcs ${params.n_pcs} --dims_use '${params.subcluster_dims}' \\
        --batch_column '${params.batch_column}' --batch_correct '${params.batch_correct}' \\
        --knn_k ${params.bcell_k} --resolution ${params.bcell_res} --n_iterations ${params.n_iterations} \\
        --top_n ${params.top_n} --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

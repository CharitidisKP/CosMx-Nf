process SUBCLUSTER {
    tag   { id }
    label 'process_highmem'
    publishDir path: { "${params.outdir}/figures/05_subclustering" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/05_subclustering" }, mode: 'copy', pattern: "*.csv"
    // One object, one folder. Each run appends its sub_<res> column to the previous one
    publishDir path: { "${params.outdir}/objects/subcluster" }, mode: 'copy', pattern: "giotto"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), emit: obj
        path "subcluster_column.txt", emit: col
        path "*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    subcluster.R --input '${obj}' --basis '${params.subcluster_basis}' \\
        --clusters '${params.subcluster_clusters}' --resolution ${params.subcluster_res} \\
        --knn_k ${params.subcluster_knn_k} --n_iterations ${params.n_iterations} \\
        --n_pcs ${params.n_pcs} --dims_use '${params.subcluster_dims}' \\
        --n_hvgs ${params.n_hvgs} --markers '${params.markers}' \\
        --batch_column '${params.batch_column}' --batch_correct '${params.batch_correct}' \\
        --min_cells ${params.subcluster_min_cells} --seed ${params.seed} \\
        --outdir . --python '${params.python_path}'
    """
}

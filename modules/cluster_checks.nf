// Checks on discover's choices: resolutions and stability by default, PCs and genes on request
process CLUSTER_CHECKS {
    tag   { id }
    label 'process_checks'
    publishDir path: "${params.outdir}/tables/03_clustering/checks", mode: 'copy',
               pattern: "{resolution_checks.csv,clustree_nodes.csv,stability_clusters.csv,picked_resolution.txt}"
    publishDir path: "${params.outdir}/figures/03_clustering/checks", mode: 'copy',
               pattern: "{clustree,resolution_checks}.png"
    publishDir path: "${params.outdir}/tables/02_integration/checks", mode: 'copy',
               pattern: "{dims,genes}_*.{csv,txt}"
    publishDir path: "${params.outdir}/figures/02_integration/checks", mode: 'copy',
               pattern: "dims_*.png"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), val(red)

    output:
        path "*.csv", emit: tables
        path "*.png"
        path "*.txt", optional: true
        path "picked_resolution.txt", emit: picked, optional: true

    script:
    def hvg = params.hvg_file ? file(params.hvg_file).toString() : ''
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} BLIS_NUM_THREADS=${task.cpus}

    cluster_checks.R --input '${obj}' --reduction '${red}' \\
        --leiden_res_ladder '${params.leiden_res_ladder}' --chosen_resolutions '${params.chosen_resolutions}' \\
        --dims_use '${params.dims_use}' --knn_k ${params.knn_k} --n_iterations ${params.n_iterations} \\
        --markers '${params.markers}' --n_hvgs ${params.n_hvgs} --hvg_method '${params.hvg_method}' \\
        --hvg_file '${hvg}' --batch_column '${params.batch_column}' --cells ${params.checks_cells} \\
        --stability ${params.checks_stability} --dims ${params.checks_dims} --genes ${params.checks_genes} \\
        --seed ${params.seed} --cores ${task.cpus} --outdir . --python '${params.python_path}'
    """
}

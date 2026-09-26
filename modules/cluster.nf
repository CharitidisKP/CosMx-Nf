process CLUSTER {
    tag   { id }
    label 'process_highmem'
    publishDir path: { id == 'merged' ? "${params.outdir}/figures/03_clustering"
                                      : "${params.outdir}/solo/${id}/figures/03_clustering" },
               mode: 'copy', pattern: "*.png"
    publishDir path: { id == 'merged' ? "${params.outdir}/tables/03_clustering"
                                      : "${params.outdir}/solo/${id}/tables/03_clustering" },
               mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), path(red)

    output:
        tuple val(id), path("giotto"), path("reduction.txt"), emit: obj
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    cluster.R --input '${obj}' --reduction_file '${red}' \\
        --knn_k ${params.knn_k} --n_iterations ${params.n_iterations} \\
        --leiden_res_ladder '${params.leiden_res_ladder}' --umap_resolutions '${params.chosen_resolutions}' \\
        --dims_use '${params.dims_use}' \\
        --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

process INTEGRATE {
    tag   { id }
    label 'process_highmem'
    publishDir path: { id == 'merged' ? "${params.outdir}/figures/02_integration"
                                      : "${params.outdir}/solo/${id}/figures/02_integration" },
               mode: 'copy', pattern: "*.png"
    publishDir path: { id == 'merged' ? "${params.outdir}/tables/02_integration"
                                      : "${params.outdir}/solo/${id}/tables/02_integration" },
               mode: 'copy', pattern: "*.{csv,txt}"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), path("reduction.txt"), emit: obj
        path "*.csv", emit: tables
        path "*.png"
        path "*.txt"

    script:
    """
    integrate.R --input '${obj}' --markers '${params.markers}' --scalefactor '${params.scalefactor}' \\
        --n_hvgs ${params.n_hvgs} --n_pcs ${params.n_pcs} --batch_column '${params.batch_column}' \\
        --batch_correct '${params.batch_correct}' \\
        --dims_use '${params.dims_use}' --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

process BANKSY {
    tag   'banksy'
    label 'process_highmem'
    publishDir path: { "${params.outdir}/figures/07_spatial/banksy" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/07_spatial/banksy" }, mode: 'copy', pattern: "*.csv"
    // Carries the BANKSY domains; --stage cci and --stage de read this one
    publishDir path: { "${params.outdir}/objects/spatial" }, mode: 'copy', pattern: "giotto"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    banksy.R --input '${obj}' --markers '${params.markers}' --n_hvgs ${params.n_hvgs} \\
        --k_geom '${params.banksy_k_geom}' --lambda_ladder '${params.banksy_lambda_ladder}' \\
        --use_agf ${params.banksy_use_agf} --n_pcs ${params.banksy_pcs} --k_nn ${params.banksy_knn} \\
        --banksy_res ${params.banksy_res} --batch_column '${params.batch_column}' \\
        --celltype_column '${params.celltype_column}' --seed ${params.seed} \\
        --outdir . --python '${params.python_path}'
    """
}

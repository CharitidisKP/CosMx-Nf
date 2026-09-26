process CCI_MISTY {
    tag   'misty'
    label 'process_cci'
    publishDir path: { "${params.outdir}/figures/08_cci/misty" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/08_cci/misty" }, mode: 'copy', pattern: "misty_*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        path "misty_*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    cci_misty.R --input '${obj}' --markers '${params.markers}' \\
        --juxta_thr ${params.misty_juxta_thr} --para_ls '${params.misty_para_ls}' \\
        --min_cells ${params.misty_min_cells} --cores ${task.cpus} --seed ${params.seed} \\
        --outdir . --python '${params.python_path}'
    """
}

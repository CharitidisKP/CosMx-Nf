// One task per sample, so every sample starts with fresh workers and -resume keeps finished ones
process CCI_MISTY {
    tag   "misty:${sample}"
    label 'process_cci'

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), val(sample)

    output:
        path "misty_*", type: 'dir', emit: results, optional: true
        path "misty_*_key.csv",      emit: keys,    optional: true

    script:
    def reuse = params.misty_reuse_dir ? "--reuse_dir '${file(params.misty_reuse_dir)}'" : ""
    """
    cci_misty.R --input '${obj}' --sample '${sample}' --markers '${params.markers}' \\
        --juxta_thr ${params.misty_juxta_thr} --para_ls '${params.misty_para_ls}' \\
        --min_cells ${params.misty_min_cells} --cores ${task.cpus} --seed ${params.seed} \\
        ${reuse} --outdir . --python '${params.python_path}'
    """
}

process CCI_MISTY_COLLECT {
    tag   'misty'
    label 'process_light'
    publishDir path: { "${params.outdir}/figures/08_cci/misty" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/08_cci/misty" }, mode: 'copy', pattern: "misty_*.csv"

    input:
        path results

    output:
        path "misty_{improvements,contributions,importances}.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    cci_misty_collect.R --input . --outdir .
    """
}

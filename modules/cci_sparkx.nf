process CCI_SPARKX {
    tag   'sparkx'
    label 'process_cci'
    publishDir path: { "${params.outdir}/figures/08_cci/sparkx" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/08_cci/sparkx" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        path "sparkx_svgs.csv", emit: svgs
        path "*.png", optional: true

    script:
    """
    cci_sparkx.R --input '${obj}' --celltype_column '${params.celltype_column}' \\
        --fdr ${params.sparkx_fdr} --cores ${task.cpus} \\
        --outdir . --python '${params.python_path}'
    """
}

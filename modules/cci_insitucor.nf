process CCI_INSITUCOR {
    tag   'insitucor'
    label 'process_cci'
    publishDir path: { "${params.outdir}/tables/08_cci/insitucor" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        path "*.csv", emit: tables

    script:
    """
    cci_insitucor.R --input '${obj}' --celltype_column '${params.celltype_column}' \\
        --k ${params.insitucor_k} --outdir . --python '${params.python_path}'
    """
}

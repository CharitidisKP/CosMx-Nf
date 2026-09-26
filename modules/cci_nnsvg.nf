process CCI_NNSVG {
    tag   'nnsvg'
    label 'process_cci'
    publishDir path: { "${params.outdir}/tables/08_cci/nnsvg" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path svg_genes

    output:
        path "nnsvg_spatially_variable_genes.csv", emit: tables

    script:
    def genes_arg = params.nnsvg_use_sparkx ? "--svg_genes '${svg_genes}'" : ""
    """
    cci_nnsvg.R --input '${obj}' ${genes_arg} --max_genes ${params.nnsvg_max_genes} \\
        --cores ${task.cpus} --outdir . --python '${params.python_path}'
    """
}

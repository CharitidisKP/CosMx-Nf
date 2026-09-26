process CCI_LIANA {
    tag   'liana'
    label 'process_cci'
    publishDir path: { "${params.outdir}/figures/08_cci/liana" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/08_cci/liana" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        path "*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    cci_liana.R --input '${obj}' --celltype_column '${params.celltype_column}' \\
        --exclude '${params.untyped_celltypes}' \\
        --min_cells ${params.liana_min_cells} --cores ${task.cpus} --seed ${params.seed} \\
        --outdir . --python '${params.python_path}'
    """
}

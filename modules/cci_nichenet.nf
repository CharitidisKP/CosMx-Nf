process CCI_NICHENET {
    tag   'nichenet'
    label 'process_cci'
    publishDir path: { "${params.outdir}/figures/08_cci/nichenet" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/08_cci/nichenet" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path samplesheet

    output:
        path "*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    cci_nichenet.R --input '${obj}' --model_dir '${params.nichenet_model_dir}' \\
        --samplesheet '${samplesheet}' --celltype_column '${params.celltype_column}' \\
        --exclude '${params.untyped_celltypes}' \\
        --receivers '${params.nichenet_receivers}' \\
        --compare_by '${params.compare_by}' --compare_baseline '${params.compare_baseline}' \\
        --compare_within '${params.compare_within}' --pair_by '${params.pair_by}' \\
        --top_ligands ${params.nichenet_top_ligands} --n_targets ${params.nichenet_n_targets} \\
        --outdir . --python '${params.python_path}'
    """
}

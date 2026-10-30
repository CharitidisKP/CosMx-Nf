process SMIDE_OVERLAP {
    tag   "${sample}"
    label 'process_smide'
    publishDir path: { "${params.outdir}/tables/09_de/smide/${sample}" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), val(sample)
        val cell_types

    output:
        path "*_smide_niche_overlap.csv", emit: ovl

    script:
    """
    smide_overlap.R --input '${obj}' --sample '${sample}' --cell_types '${cell_types}' \\
        --celltype_column '${params.celltype_column}' --groupvar '${params.smide_groupvar}' \\
        --radius ${params.smide_radius} --outdir . --python '${params.python_path}'
    """
}

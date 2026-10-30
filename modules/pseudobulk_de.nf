process PSEUDOBULK_DE {
    tag   'pseudobulk_de'
    label 'process_diag'
    publishDir path: { "${params.outdir}/figures/09_de" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/09_de/pseudobulk" }, mode: 'copy', pattern: "pb_*.{csv,csv.gz}"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path samplesheet

    output:
        path "pb_{summary,de}.csv", emit: tables
        path "pb_{design,per_subject_lfc}.csv"
        path "pb_counts.csv.gz"
        path "*.png", optional: true

    script:
    """
    pseudobulk_de.R --input '${obj}' --samplesheet '${samplesheet}' \\
        --compare_by '${params.compare_by}' --compare_baseline '${params.compare_baseline}' \\
        --compare_within '${params.compare_within}' --pair_by '${params.pair_by}' \\
        --celltype_column '${params.celltype_column}' --exclude '${params.untyped_celltypes}' \\
        --min_cells ${params.pb_min_cells} --fov_min_cells ${params.pb_fov_min_cells} \\
        --fdr ${params.pb_fdr} --interaction ${params.pb_interaction} \\
        --radius ${params.smide_radius} --overlap_flag ${params.pb_overlap_flag} \\
        --outdir . --python '${params.python_path}'
    """
}

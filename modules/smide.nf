process SMIDE {
    tag   "${sample}:${cell_type}"
    label 'process_smide'
    publishDir path: { "${params.outdir}/tables/09_de/smide/${sample}" }, mode: 'copy', pattern: "*.{csv,txt}"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), val(sample), val(cell_type)

    output:
        path "*_smide_de_results.csv",   emit: de,      optional: true
        path "*_smide_overlap_ratio.csv", emit: overlap, optional: true
        path "*_smide_skipped.csv",       emit: skipped
        path "*_smide_failed_genes.txt",  emit: failed,  optional: true

    script:
    def fmla = params.smide_formula ? "--formula '${params.smide_formula}'" : ""
    """
    export OMP_NUM_THREADS=1 OPENBLAS_NUM_THREADS=1 MKL_NUM_THREADS=1 BLIS_NUM_THREADS=1

    smide.R --input '${obj}' --sample '${sample}' --cell_type '${cell_type}' \\
        --celltype_column '${params.celltype_column}' --groupvar '${params.smide_groupvar}' \\
        --radius ${params.smide_radius} --overlap_threshold ${params.smide_overlap_threshold} \\
        --min_detection ${params.smide_min_detection} \\
        --min_cells_per_niche ${params.smide_min_cells_per_niche} \\
        --family ${params.smide_family} --spatial_model '${params.smide_spatial_model}' \\
        --k_prop_n ${params.smide_k_prop_n} --comparisons '${params.smide_comparisons}' \\
        ${fmla} --cores ${task.cpus} --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

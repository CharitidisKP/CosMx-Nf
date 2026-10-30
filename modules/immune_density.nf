process IMMUNE_DENSITY {
    tag   'immune_density'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/figures/07_spatial/density" }, mode: 'copy', pattern: "*.{png,pdf,jpg}"
    publishDir path: { "${params.outdir}/tables/07_spatial/density" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path samplesheet

    output:
        path "immune_density.csv", emit: tables
        path "*.{png,pdf,jpg}"

    script:
    """
    immune_density.R --input '${obj}' --samplesheet '${samplesheet}' \\
        --compare_by '${params.compare_by}' --compare_baseline '${params.compare_baseline}' \\
        --compare_within '${params.compare_within}' --pair_by '${params.pair_by}' \\
        --celltype_column '${params.celltype_column}' --immune_types '${params.immune_types}' \\
        --um_per_px ${params.um_per_px} --bin_um ${params.density_bin_um} \\
        --label_colours '${params.label_colours}' \\
        --fig_format '${params.fig_format}' --outdir . --python '${params.python_path}'
    """
}

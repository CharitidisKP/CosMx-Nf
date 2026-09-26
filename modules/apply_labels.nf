process APPLY_LABELS {
    tag   'apply_labels'
    label 'process_light'
    publishDir path: { "${params.outdir}/figures/06_annotation" }, mode: 'copy', pattern: "*.png"
    publishDir path: { "${params.outdir}/tables/06_annotation" }, mode: 'copy', pattern: "*.csv"
    // The labelled object is what --stage spatial picks up by default
    publishDir path: { "${params.outdir}/objects/annotate" }, mode: 'copy', pattern: "giotto"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), path(labels)

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    apply_labels.R --input '${obj}' --labels '${labels}' --cluster_basis '${params.cluster_basis}' \\
        --immune_types '${params.immune_types}' --markers '${params.markers}' \\
        --outdir . --python '${params.python_path}'
    """
}

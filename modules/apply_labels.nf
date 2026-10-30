process APPLY_LABELS {
    tag   'apply_labels'
    label 'process_light'
    publishDir path: { "${params.outdir}/figures/06_annotation" }, mode: 'copy', pattern: "*.{png,pdf,jpg}"
    publishDir path: { "${params.outdir}/tables/06_annotation" }, mode: 'copy', pattern: "*.csv"
    // The labelled object is what --stage spatial picks up by default, with the label state it
    // was built from beside it
    publishDir path: { "${params.outdir}/objects/annotate" }, mode: 'copy', pattern: "{giotto,label_state.txt}"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), path(labels, stageAs: 'cluster_labels.csv'), path(overrides, stageAs: 'label_overrides.csv'), val(st)

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "label_state.txt"
        path "*.{png,pdf,jpg}"

    script:
    """
    apply_labels.R --input '${obj}' --labels '${labels}' --overrides '${overrides}' \\
        --cluster_basis '${st.basis}' --label_state '${st.digest}' \\
        --labels_from '${st.labels_from}' --overrides_from '${st.overrides_from}' \\
        --immune_types '${params.immune_types}' --untyped '${params.untyped_celltypes}' \\
        --markers '${params.markers}' --label_colours '${params.label_colours}' \\
        --fig_format '${params.fig_format}' --outdir . --python '${params.python_path}'
    """
}

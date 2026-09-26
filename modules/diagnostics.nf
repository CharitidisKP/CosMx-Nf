process DIAGNOSTICS {
    tag   { "${id}:${col}" }
    label 'process_diag'
    // Each column lands with the step that made it: clusters, reference calls or subclusters
    publishDir path: {
        def step = col.startsWith('sub_') ? '05_subclustering'
                 : (col.startsWith('insitutype') || col.startsWith('ht_')) ? '04_reference_typing'
                 : '03_clustering'
        id == 'merged' ? "${params.outdir}/tables/${step}/markers"
                       : "${params.outdir}/solo/${id}/tables/${step}/markers"
    }, mode: 'copy', pattern: "*.csv"
    publishDir path: {
        def step = col.startsWith('sub_') ? '05_subclustering'
                 : (col.startsWith('insitutype') || col.startsWith('ht_')) ? '04_reference_typing'
                 : '03_clustering'
        id == 'merged' ? "${params.outdir}/figures/${step}"
                       : "${params.outdir}/solo/${id}/figures/${step}"
    }, mode: 'copy', pattern: "*.png"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), val(col)

    output:
        path "*.csv", emit: tables
        path "*.png", optional: true

    script:
    """
    diagnostics.R --input '${obj}' --cluster_column '${col}' --markers '${params.markers}' \\
        --top_n ${params.top_n} --min_detection ${params.marker_min_detection} \\
        --min_panel_genes ${params.marker_min_panel_genes} --min_score ${params.marker_min_score} \\
        --outdir . --python '${params.python_path}'
    """
}

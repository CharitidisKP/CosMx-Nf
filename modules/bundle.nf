// Writes the Shiny app's bundle, local_export/, from the labelled object
process BUNDLE {
    tag   { id }
    label 'process_load'
    publishDir path: { "${params.outdir}/local_export" }, mode: 'copy', overwrite: true

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        // The cells table also goes to REPORT, which draws the UMAPs from it
        path "cells_umap_clusters.rds", emit: cells
        path "*.rds"
        path "export_app_bundle.txt"

    script:
    """
    export_bundle.R --input '${obj}' --outdir . --log export_app_bundle.txt \\
        --genes all --poly_compress xz --markers '${params.markers}' \\
        --tables '${params.outdir}/tables' --python '${params.python_path}'
    """
}

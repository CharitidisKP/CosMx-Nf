// Scores every labelling against FOVs, cell types and space, and draws the maps
process BANKSY_SWEEP_COLLECT {
    tag   'collect'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/tables/07_spatial/banksy_sweep" },  mode: 'copy', pattern: "*.csv"
    publishDir path: { "${params.outdir}/figures/07_spatial/banksy_sweep" }, mode: 'copy', pattern: "*.png"

    input:
        path(labels, stageAs: 'labels/*')
        path(times,  stageAs: 'times/*')
        path(metas,  stageAs: 'meta/*')

    output:
        path "*.csv", emit: tables
        path "*.png"

    script:
    """
    banksy_sweep_collect.R --labels labels --times times --meta meta \\
        --map_res '${params.banksy_sweep_map_res}' --untyped '${params.untyped_celltypes}' \\
        --focal '${params.banksy_sweep_focal}' --outdir .
    """
}

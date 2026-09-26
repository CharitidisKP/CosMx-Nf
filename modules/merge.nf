process MERGE {
    tag   'merge'
    label 'process_highmem'
    publishDir path: { "${params.outdir}/tables/01_qc" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(sids), path(dirs, stageAs: 'input_giotto_*')

    output:
        tuple val('merged'), path("giotto"), emit: obj
        path "*.csv", emit: tables

    script:
    """
    merge.R --sample_ids '${sids.join(',')}' --x_padding ${params.x_padding} --outdir . --python '${params.python_path}'
    """
}

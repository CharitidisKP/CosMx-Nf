process LOAD {
    tag   { meta.id }
    label 'process_load'
    publishDir path: { "${params.outdir}/tables/01_qc" }, mode: 'copy', pattern: "*.csv"

    input:
        tuple val(meta), path(data_dir)

    output:
        tuple val(meta), path("giotto"), emit: obj
        path "*.csv", emit: tables

    script:
    """
    load.R --sample_id '${meta.id}' --file_prefix '${meta.prefix}' --data_dir '${data_dir}' \\
        --patient_id '${meta.patient}' --treatment '${meta.treatment}' --timepoint '${meta.timepoint}' \\
        --slide_id '${meta.slide}' --batch '${meta.batch}' \\
        --subset_fovs '${meta.subset_fovs}' --outdir . --python '${params.python_path}'
    """
}

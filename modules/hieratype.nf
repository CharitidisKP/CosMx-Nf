process HIERATYPE {
    tag   { id }
    label 'process_highmem'
    publishDir path: { id == 'merged' ? "${params.outdir}/tables/04_reference_typing/hieratype"
                                      : "${params.outdir}/solo/${id}/tables/04_reference_typing/hieratype" },
               mode: 'copy', pattern: "*.csv"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto'), path(red)

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables

    script:
    """
    hieratype.R --input '${obj}' --reduction_file '${red}' \\
        --dims_use '${params.dims_use}' --knn_k ${params.hieratype_knn_k} \\
        --batch_column '${params.hieratype_batch_column}' \\
        --call_top ${params.hieratype_call_top} --call_second ${params.hieratype_call_second} \\
        --seed ${params.seed} --outdir . --python '${params.python_path}'
    """
}

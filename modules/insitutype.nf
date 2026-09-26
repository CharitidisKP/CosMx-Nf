process INSITUTYPE {
    tag   { id }
    label 'process_insitutype'
    publishDir path: { id == 'merged' ? "${params.outdir}/tables/04_reference_typing/insitutype_unsup"
                                      : "${params.outdir}/solo/${id}/tables/04_reference_typing/insitutype_unsup" },
               mode: 'copy', pattern: "*.csv"
    // With supervised profiles, ATTACH_SUPERVISED publishes the discover object instead
    publishDir path: { "${params.outdir}/objects" }, mode: 'copy', pattern: "giotto",
               saveAs: { _fn -> id != 'merged' ? "solo/${id}/giotto"
                                              : (params.supervised_profiles ? null : "discover/giotto") }

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables

    script:
    """
    insitutype.R --input '${obj}' --auto_lo ${params.insitutype_auto_lo} --auto_hi ${params.insitutype_auto_hi} \\
        --n_starts ${params.insitutype_n_starts} --seed ${params.seed} \\
        --cohort ${params.insitutype_cohort} --cohort_vars '${params.cohort_vars}' \\
        --cohort_column '${params.insitutype_cohort_column ?: ""}' \\
        --nb_k ${params.insitutype_nb_k} --nb_pcs ${params.insitutype_nb_pcs} \\
        --reference '${params.insitutype_reference ?: ""}' \\
        --semi_basis '${params.insitutype_semi_basis ?: ""}' --semi_n '${params.insitutype_semi_n ?: ""}' \\
        --outdir . --python '${params.python_path}'
    """
}

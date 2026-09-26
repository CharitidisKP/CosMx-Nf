// One task per reference profile, emitting per-cell CSVs; ATTACH_SUPERVISED folds them onto
// the Giotto object in a single write.
process SUPERVISED {
    tag { "${id}:${prof.name}" }
    label 'process_insitutype'
    publishDir path: {
        id == 'merged'
            ? "${params.outdir}/tables/04_reference_typing/supervised"
            : "${params.outdir}/solo/${id}/tables/04_reference_typing/supervised"
    }, mode: 'copy', pattern: '*.csv'

    input:
    tuple val(id), path(obj, stageAs: 'input_giotto'), val(prof)

    output:
    tuple val(id), path("insitutype_sup_${prof.name}_cells.csv"), emit: cells
    path "*.csv", emit: tables

    script:
    """
      insitutype_supervised.R --input '${obj}' --reference '${prof.reference}' \\
          --name '${prof.name}' --refine ${params.insitutype_refine} \\
          --conf_threshold ${params.insitutype_conf_threshold} \\
          --min_gene_overlap ${params.supervised_min_gene_overlap} \\
          --cohort ${params.insitutype_cohort} --cohort_vars '${params.cohort_vars}' \\
          --seed ${params.seed} --outdir . --python '${params.python_path}'
      """
}

process ATTACH_SUPERVISED {
    tag { id }
    label 'process_load'
    publishDir path: {
        id == 'merged'
            ? "${params.outdir}/objects/discover"
            : "${params.outdir}/objects/solo/${id}"
    }, mode: 'copy', pattern: 'giotto'

    input:
    tuple val(id), path(obj, stageAs: 'input_giotto'), path(csvs)

    output:
    tuple val(id), path("giotto"), emit: obj

    script:
    """
    attach_supervised.R --input '${obj}' --outdir . --python '${params.python_path}'
    """
}

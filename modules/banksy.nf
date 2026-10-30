// Load the object once: pooled variable genes, and per sample its counts and cells.
// samples: the samples to write, blank for every one
process BANKSY_PREP {
    tag   'prep'
    label 'process_spatial'

    input:
        path(obj, stageAs: 'input_giotto')
        val samples

    output:
        path "hvg_pooled.txt", emit: hvg
        path "*_raw.rds",      emit: raw
        path "*_meta.csv.gz",  emit: meta

    script:
    """
    banksy_prep.R --input '${obj}' --samples '${samples}' \\
        --hvgs ${params.banksy_hvgs} --celltype_column '${params.celltype_column}' \\
        --outdir . --python '${params.python_path}'
    """
}

// One sample through BANKSY with one settings set, v, a row of banksy_variants. The spatial
// stage runs it with the row banksy_settings names, at the pipeline's lambdas and resolution.
// The sweep runs it as BANKSY_SWEEP with every row; conf/base.config sets what differs there
process BANKSY_SAMPLE {
    tag   "${sample}"
    label 'process_sweep'

    input:
        tuple val(sample), path(raw), path(meta), val(v)
        path(hvg)

    output:
        path "*_labels.csv.gz", emit: labels
        path "*_genes.csv",     emit: genes
        path "*_time.csv",      emit: time

    script:
    """
    export OMP_NUM_THREADS=${task.cpus} OPENBLAS_NUM_THREADS=${task.cpus} MKL_NUM_THREADS=${task.cpus} BLIS_NUM_THREADS=${task.cpus}

    banksy_sample.R --raw '${raw}' --meta '${meta}' --hvg '${hvg}' --sample '${sample}' \\
        --variant '${v.variant}' --genes '${v.genes}' --values '${v.values}' \\
        --k_geom '${v.k_geom}' --agf ${v.agf} --npcs ${v.npcs} --group '${v.group}' \\
        --lambdas '${v.lambdas}' --resolutions '${v.resolutions}' --k_nn ${v.k_nn} \\
        --cluster '${v.cluster}' --bg_mult ${params.banksy_bg_mult} \\
        --untyped '${params.untyped_celltypes}' --seed ${params.seed} --outdir .
    """
}

// Each sample's niches matched into types across samples, written to the object, counted and drawn
process BANKSY {
    tag   'banksy'
    label 'process_spatial'
    publishDir path: { "${params.outdir}/figures/07_spatial/banksy" }, mode: 'copy', pattern: "*.{png,pdf,jpg}"
    publishDir path: { "${params.outdir}/tables/07_spatial/banksy" }, mode: 'copy', pattern: "*.csv"
    // Carries the BANKSY domains; --stage cci and --stage de read this one, with the label
    // state it was built from beside it
    publishDir path: { "${params.outdir}/objects/spatial" }, mode: 'copy', pattern: "{giotto,label_state.txt}"

    input:
        tuple val(id), path(obj, stageAs: 'input_giotto')
        path(labels, stageAs: 'labels/*')
        path(genes, stageAs: 'genes/*')
        val v
        val state

    output:
        tuple val(id), path("giotto"), emit: obj
        path "*.csv", emit: tables
        path "label_state.txt"
        path "*.{png,pdf,jpg}"

    script:
    """
    banksy.R --input '${obj}' --labels labels --genes genes --variant '${v.variant}' \\
        --k_geom '${v.k_geom}' --use_agf ${v.agf} --k_nn ${v.k_nn} \\
        --lambda_ladder '${v.lambdas}' --sample_res ${v.resolutions} \\
        --match_cut ${params.banksy_match_cut} --bg_mult ${params.banksy_bg_mult} \\
        --untyped '${params.untyped_celltypes}' --celltype_column '${params.celltype_column}' \\
        --um_per_px ${params.um_per_px} --focus '${params.niche_map_focus}' \\
        --label_colours '${params.label_colours}' \\
        --fig_format '${params.fig_format}' --outdir . --python '${params.python_path}'
    printf '%s\\n' '${state}' > label_state.txt
    """
}

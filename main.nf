nextflow.enable.dsl = 2

include { DISCOVERY }   from './subworkflows/discovery'
include { ANNOTATION; SPATIAL; DE; banksy_settings; banksy_key } from './subworkflows/annotation'
include { CCI }         from './subworkflows/cci'
include { DIAGNOSTICS } from './modules/diagnostics'
include { SUPERVISED; ATTACH_SUPERVISED } from './modules/insitutype_supervised'
include { SUBCLUSTER }  from './modules/subcluster'
include { REPORT }      from './modules/report'
// The sweep runs the spatial stage's BANKSY steps under its own names, with every settings set
include { BANKSY_PREP as BANKSY_SWEEP_PREP; BANKSY_SAMPLE as BANKSY_SWEEP } from './modules/banksy'
include { BANKSY_SWEEP_COLLECT } from './modules/banksy_sweep'
include { isIncluded }  from './subworkflows/samplesheet'

// Short commit of the pipeline checkout, recorded in run_report.md
def git_commit() {
    try {
        return "git -C ${projectDir} rev-parse --short HEAD".execute().text.trim() ?: "unknown"
    } catch (Exception _e) {
        return "unknown"
    }
}

workflow discover {
    def sheet_ids = file(params.input).splitCsv(header: true).findAll { isIncluded(it) }.collect { it.sample_id }
    def solo_ids  = !params.solo ? [] :
        (params.solo.toString().trim() == 'all' ? sheet_ids
                                                : params.solo.toString().split(',').collect { it.trim() }.findAll { it })
    def unknown = solo_ids - sheet_ids
    if (unknown)
        error "--solo: unknown sample_id(s): ${unknown.join(', ')}. Available: ${sheet_ids.join(', ')}"
    if (params.solo_only && !solo_ids)
        error "--solo_only requires --solo <sample_id[,sample_id,...]>"

    DISCOVERY( file(params.input), solo_ids )
    REPORT( params.stage, git_commit(), DISCOVERY.out.tables.collect().ifEmpty([]), file(params.input) )
}

// Each stage picks up where the last one published, so they chain without flags.
def stage_obj(stage) {
    params.merged_giotto ?: "${params.outdir}/objects/${stage}/giotto"
}

// The run's curated labels, beside its results: cluster_labels.<basis>.csv for each clustering
// labelled, label_overrides.csv with the per-cell edits, and cluster_basis.txt naming the
// clustering cell_type comes from. The Shiny app writes the same files
def labels_dir() {
    file(params.labels_dir ?: "${params.outdir}/../labels").normalize()
}

// What annotate applies: the clustering, its cluster labels and the per-cell edits, each read
// from labels_dir unless a param names it, and a digest of the three. A stage compares the
// digest with the one its input object was built from
def label_state() {
    def dir = labels_dir()
    def basis = params.cluster_basis ? params.cluster_basis.toString() : null
    def pick = file("${dir}/cluster_basis.txt")
    if (!basis && pick.exists()) basis = pick.text.trim()
    if (!basis) {
        def found = dir.isDirectory() ? dir.list().findAll { n -> n ==~ /cluster_labels\..+\.csv/ } : []
        def bases = found.collect { n -> n.replaceFirst(/^cluster_labels\./, '').replaceFirst(/\.csv$/, '') }
        if (bases.size() != 1)
            error "No clustering chosen for cell_type. Write its name to ${dir}/cluster_basis.txt " +
                  "(the app's Cluster labels panel does), or pass --cluster_basis. " +
                  (bases ? "Labels found for: ${bases.join(', ')}" : "${dir} holds no cluster_labels.<basis>.csv")
        basis = bases[0]
    }
    def labels = file(params.labels ?: "${dir}/cluster_labels.${basis}.csv")
    if (!labels.exists())
        error "No cluster labels for ${basis} at ${labels}. Save them in the app's Cluster labels panel, " +
              "copy the draft from tables/06_annotation/labels/, or pass --labels."
    def overrides = file(params.label_overrides ?: "${dir}/label_overrides.csv")
    if (params.label_overrides && !overrides.exists())
        error "No label overrides at ${overrides}"
    // Without edits, a file holding only the header
    if (!overrides.exists()) overrides = file("${projectDir}/assets/no_label_overrides.csv")
    def digest = "${basis}\n${labels.text}\n${overrides.text}".toString().md5()
    [basis: basis, labels: labels, overrides: overrides, digest: digest]
}

// The label state a stage's published object was built from, or null
def built_from(stage) {
    def f = file("${params.outdir}/objects/${stage}/label_state.txt")
    f.exists() ? f.text.trim() : null
}

// The clustering object annotate reads: subcluster columns live on the subcluster object,
// Leiden columns on the discover one
def clustering_obj(basis) {
    def mg = stage_obj(basis.startsWith('sub_') ? 'subcluster' : 'discover')
    if (!file(mg).exists())
        error "annotate: no Giotto object at ${mg}. Pass --merged_giotto <dir> (a saveGiotto folder)."
    mg
}

// Attach cell_type from the run's labels, and count it.
workflow annotate {
    def st = label_state()
    log.info "annotate: ${st.basis}, labels ${st.labels}, edits ${st.overrides}"
    ANNOTATION( clustering_obj(st.basis), st )
    REPORT( params.stage, git_commit(), ANNOTATION.out.tables.collect().ifEmpty([]), file(params.input) )
}

// The labelled object for a later stage. The published one when it was built from the current
// labels; otherwise annotate runs first, in this launch. --merged_giotto is taken as it is
workflow LABELLED {
    main:
        def st = label_state()
        // Not def: emit reads them
        obj = Channel.empty()
        counts = Channel.empty()
        tables = Channel.empty()
        if (params.merged_giotto || built_from('annotate') == st.digest) {
            def mg = stage_obj('annotate')
            if (!file(mg).exists())
                error "No labelled object at ${mg}. Run --stage annotate, or pass --merged_giotto."
            obj = Channel.value(tuple('merged', file(mg)))
            counts = Channel.fromPath("${params.outdir}/tables/06_annotation/cell_type_counts.csv")
        } else {
            log.info "The labels changed since annotate last ran (or it recorded none), so annotate runs first"
            ANNOTATION( clustering_obj(st.basis), st )
            obj = ANNOTATION.out.annotated
            counts = ANNOTATION.out.counts
            tables = ANNOTATION.out.tables
        }

    emit:
        obj    = obj
        counts = counts
        tables = tables
        digest = Channel.value(st.digest)
}

// The spatial object for cci and de, rebuilt in this launch when it is behind the labels or the
// BANKSY settings, with annotate first when that is behind the labels too
workflow SPATIAL_OBJECT {
    main:
        def st = label_state()
        // Not def: emit reads them
        obj = Channel.empty()
        counts = Channel.empty()
        tables = Channel.empty()
        if (params.merged_giotto || built_from('spatial') == "${st.digest}:${banksy_key(banksy_settings())}".toString()) {
            def mg = stage_obj('spatial')
            if (!file(mg).exists())
                error "No spatial object at ${mg}. Run --stage spatial, or pass --merged_giotto."
            obj = Channel.value(tuple('merged', file(mg)))
            counts = Channel.fromPath("${params.outdir}/tables/06_annotation/cell_type_counts.csv")
        } else {
            log.info "The labels or the BANKSY settings changed since the spatial stage last ran (or it recorded none), so it runs first"
            LABELLED()
            SPATIAL( LABELLED.out.obj, LABELLED.out.digest )
            obj = SPATIAL.out.spatial
            counts = LABELLED.out.counts
            tables = LABELLED.out.tables.mix(SPATIAL.out.tables)
        }

    emit:
        obj    = obj
        counts = counts
        tables = tables
}

// Delaunay network, BANKSY domains, B-lineage subclustering.
workflow spatial {
    LABELLED()
    SPATIAL( LABELLED.out.obj, LABELLED.out.digest )
    REPORT( params.stage, git_commit(),
            LABELLED.out.tables.mix(SPATIAL.out.tables).collect().ifEmpty([]), file(params.input) )
}

// Ligand-receptor and spatial-correlation methods, on the BANKSY object.
workflow cci {
    // Checked here, so a missing model file stops the stage before any task starts
    if (params.nichenet_model_dir) {
        def nn_dir   = file(params.nichenet_model_dir)
        def nn_files = nn_dir.isDirectory() ? nn_dir.list().toList() : []
        def nn_miss  = ['ligand_target_matrix', 'lr_network', 'weighted_networks']
            .findAll { p -> !nn_files.any { n -> n.startsWith(p) && n.endsWith('.rds') } }
        if (nn_miss)
            error "cci: no ${nn_miss.join(', ')} .rds in ${params.nichenet_model_dir}. " +
                  "Stage the NicheNet model files there, or set nichenet_model_dir = '' in the config to skip NicheNet."
    }
    SPATIAL_OBJECT()
    CCI( SPATIAL_OBJECT.out.obj )
    REPORT( params.stage, git_commit(),
            SPATIAL_OBJECT.out.tables.mix(CCI.out.tables).collect().ifEmpty([]), file(params.input) )
}

// smiDE, fanned out over sample x cell type, and pseudobulk DE across samples. The cell types
// come from smide_celltypes, or else from cell_type itself
workflow de {
    SPATIAL_OBJECT()
    DE( SPATIAL_OBJECT.out.obj, SPATIAL_OBJECT.out.counts )
    REPORT( params.stage, git_commit(),
            SPATIAL_OBJECT.out.tables.mix(DE.out.tables).collect().ifEmpty([]), file(params.input) )
}

// BANKSY per biopsy over every settings set in banksy_variants, on the spatial object, by the
// same script the spatial stage runs. Publishes beside the spatial stage's BANKSY outputs and
// changes no object.
workflow banksy_sweep {
    SPATIAL_OBJECT()
    BANKSY_SWEEP_PREP( SPATIAL_OBJECT.out.obj.map { _id, obj -> obj }, params.banksy_sweep_samples )
    def ch_raw = BANKSY_SWEEP_PREP.out.raw.flatten()
        .map { f -> tuple(f.name.replaceFirst('_raw\\.rds$', ''), f) }
    def ch_meta = BANKSY_SWEEP_PREP.out.meta.flatten()
        .map { f -> tuple(f.name.replaceFirst('_meta\\.csv\\.gz$', ''), f) }
    def ch_variants = Channel.fromPath(params.banksy_variants).splitCsv(header: true)
    BANKSY_SWEEP( ch_raw.join(ch_meta).combine(ch_variants), BANKSY_SWEEP_PREP.out.hvg.first() )
    BANKSY_SWEEP_COLLECT( BANKSY_SWEEP.out.labels.collect(), BANKSY_SWEEP.out.time.collect(),
                          BANKSY_SWEEP_PREP.out.meta.collect() )
}

// Optional pass between discover and annotate. Publishes to objects/subcluster, leaving the
// discover object alone. Annotate picks it up for any sub_* basis.
workflow subcluster {
    if (!params.subcluster_basis)
        error "subcluster requires --subcluster_basis <cluster column, e.g. leiden_clus_res0.3>"
    // Build on the previous subcluster object when there is one, so a second resolution
    // appends its column instead of starting from discover again.
    def prev = file("${params.outdir}/objects/subcluster/giotto")
    def mg = params.merged_giotto ?: (prev.exists() ? prev : "${params.outdir}/objects/discover/giotto")
    if (!file(mg).exists())
        error "subcluster: no Giotto object at ${mg}. Pass --merged_giotto <dir> (a saveGiotto folder)."

    def which = params.subcluster_clusters ?: 'every cluster'
    log.info "subcluster: ${mg} -> ${params.subcluster_basis}, ${which}, resolution ${params.subcluster_res}"

    SUBCLUSTER( Channel.value(tuple('merged', file(mg))) )

    // Same triage discover produces, so subclusters arrive with a label template
    DIAGNOSTICS( SUBCLUSTER.out.obj
                   .combine(SUBCLUSTER.out.col.map { f -> f.text.trim() }) )
    REPORT( params.stage, git_commit(),
            SUBCLUSTER.out.tables.mix(DIAGNOSTICS.out.tables).collect().ifEmpty([]), file(params.input) )
}

// Runs DIAGNOSTICS against an object that already exists, without re-running
// INTEGRATE/CLUSTER/INSITUTYPE and without republishing an object.
workflow diagnose {
    def cols = params.diag_columns
        ? params.diag_columns.toString().split(',').collect { it.trim() }.findAll { it }
        : params.chosen_resolutions.toString().split(',').collect { "leiden_clus_res${it.trim()}" } + ['insitutype_unsup']
    if (!cols) error "diagnose: --diag_columns resolved to nothing."
    // Subcluster columns live on the subcluster object, which also carries every discover column
    def mg = stage_obj(cols.any { it.startsWith('sub_') } ? 'subcluster' : 'discover')
    if (!file(mg).exists())
        error "diagnose: no Giotto object at ${mg}. Pass --merged_giotto <dir> (it must be a saveGiotto folder)."
    log.info "diagnose: ${mg} -> columns ${cols.join(', ')}"

    DIAGNOSTICS( Channel.value(tuple('merged', file(mg))).combine(Channel.fromList(cols)) )
    REPORT( params.stage, git_commit(), DIAGNOSTICS.out.tables.collect().ifEmpty([]), file(params.input) )
}

// Runs the supervised passes and their diagnostics against an object that already exists,
// without re-running LOAD/QC/MERGE/INTEGRATE/CLUSTER/INSITUTYPE.
workflow supervised {
    def mg = stage_obj('discover')
    if (!file(mg).exists())
        error "supervised: no Giotto object at ${mg}. Pass --merged_giotto <dir> (it must be a saveGiotto folder)."

    // Profiles are named in the config and staged per machine under profile_dir.
    def prof_names = (params.supervised_profiles ?: "").toString()
        .split(',').collect { n -> n.trim() }.findAll { n -> n }
    if (!prof_names) error "supervised: params.supervised_profiles is empty; nothing to run."
    def profs = prof_names.collect { n -> [
        name      : n,
        reference : ['rds', 'csv']
            .collect { ext -> file("${params.profile_dir}/${n}.${ext}") }
            .find { f -> f.exists() }
    ] }
    def missing = profs.findAll { p -> !p.reference }
    if (missing)
        error "supervised: no .rds or .csv in ${params.profile_dir} for:\n" +
              missing.collect { p -> "  ${p.name}" }.join('\n')

    def sup_cols = profs.collect { p -> params.insitutype_refine
        ? ["insitutype_sup_${p.name}", "insitutype_sup_${p.name}_refined"]
        : ["insitutype_sup_${p.name}"] }.flatten()
    def cols = params.diag_columns
        ? params.diag_columns.toString().split(',').collect { it.trim() }.findAll { it }
        : sup_cols
    log.info "supervised: ${mg} -> profiles ${profs.collect{ it.name }.join(', ')}; diagnostics on ${cols.join(', ')}"

    SUPERVISED( Channel.fromList(profs).map { p -> tuple('merged', file(mg), p) } )
    ATTACH_SUPERVISED( SUPERVISED.out.cells.map { id, csv -> csv }.collect()
                         .map { csvs -> tuple('merged', file(mg), csvs) } )
    def ch_tables = SUPERVISED.out.tables
    if (cols) {
        DIAGNOSTICS( ATTACH_SUPERVISED.out.obj.combine(Channel.fromList(cols)) )
        ch_tables = ch_tables.mix(DIAGNOSTICS.out.tables)
    }
    REPORT( params.stage, git_commit(), ch_tables.collect().ifEmpty([]), file(params.input) )
}

// Rebuilds run_report.md from what is already published, e.g. after a failed stage
workflow rebuild_report {
    REPORT( params.stage, git_commit(), Channel.value([]), file(params.input) )
}

workflow {
    // A mistyped stage would otherwise fall through to discover and replace its outputs
    def stages = ['discover', 'subcluster', 'diagnose', 'supervised', 'annotate', 'spatial', 'cci', 'de',
                  'banksy_sweep', 'report']
    if (!(params.stage in stages)) error "--stage must be one of: ${stages.join(', ')}"

    // Every launch leaves its parameters and a line in the stage log beside the Nextflow reports
    def info = file("${params.outdir}/pipeline_info")
    info.mkdirs()
    file("${info}/${params.stage}_params.json").text =
        groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(params))
    def stage_log = file("${info}/stage_log.tsv")
    if (!stage_log.exists()) stage_log.text = "stage\tstarted\tcommit\tcommand\n"
    stage_log.append("${params.stage}\t${new Date().format('yyyy-MM-dd HH:mm')}\t${git_commit()}\t${workflow.commandLine}\n")

    if      (params.stage == 'annotate')   annotate()
    else if (params.stage == 'spatial')    spatial()
    else if (params.stage == 'cci')        cci()
    else if (params.stage == 'de')         de()
    else if (params.stage == 'subcluster') subcluster()
    else if (params.stage == 'diagnose')   diagnose()
    else if (params.stage == 'supervised') supervised()
    else if (params.stage == 'banksy_sweep') banksy_sweep()
    else if (params.stage == 'report')     rebuild_report()
    else                                   discover()
}

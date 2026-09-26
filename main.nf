nextflow.enable.dsl = 2

include { DISCOVERY }   from './subworkflows/discovery'
include { ANNOTATION; SPATIAL; DE } from './subworkflows/annotation'
include { CCI }         from './subworkflows/cci'
include { DIAGNOSTICS } from './modules/diagnostics'
include { SUPERVISED; ATTACH_SUPERVISED } from './modules/insitutype_supervised'
include { SUBCLUSTER }  from './modules/subcluster'
include { REPORT }      from './modules/report'
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

// Attach cell_type from the curated labels CSV, and count it.
workflow annotate {
    if (!params.labels)        error "annotate requires --labels <cluster_labels.csv>"
    if (!params.cluster_basis) error "annotate requires --cluster_basis <cluster column, e.g. sub_res0.3>"
    // Subcluster columns live on the subcluster object, Leiden columns on the discover one
    def mg = stage_obj(params.cluster_basis.toString().startsWith('sub_') ? 'subcluster' : 'discover')
    if (!file(mg).exists())
        error "annotate: no Giotto object at ${mg}. Pass --merged_giotto <dir> (a saveGiotto folder)."
    ANNOTATION( mg, params.labels )
    REPORT( params.stage, git_commit(), ANNOTATION.out.tables.collect().ifEmpty([]), file(params.input) )
}

// Delaunay network, BANKSY domains, B-lineage subclustering.
workflow spatial {
    def mg = stage_obj('annotate')
    if (!file(mg).exists())
        error "spatial: no labelled object at ${mg}. Run --stage annotate first, or pass --merged_giotto."
    SPATIAL( mg )
    REPORT( params.stage, git_commit(), SPATIAL.out.tables.collect().ifEmpty([]), file(params.input) )
}

// Ligand-receptor and spatial-correlation methods, on the BANKSY object.
workflow cci {
    def mg = stage_obj('spatial')
    if (!file(mg).exists())
        error "cci: no spatial object at ${mg}. Run --stage spatial first, or pass --merged_giotto."
    CCI( Channel.value(tuple('merged', file(mg))) )
    REPORT( params.stage, git_commit(), CCI.out.tables.collect().ifEmpty([]), file(params.input) )
}

// smiDE, fanned out over sample x cell type.
workflow de {
    if (!params.labels) error "de requires --labels <cluster_labels.csv> for the cell-type list"
    def mg = stage_obj('spatial')
    if (!file(mg).exists())
        error "de: no spatial object at ${mg}. Run --stage spatial first, or pass --merged_giotto."
    DE( mg, params.labels )
    REPORT( params.stage, git_commit(), DE.out.tables.collect().ifEmpty([]), file(params.input) )
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
    def stages = ['discover', 'subcluster', 'diagnose', 'supervised', 'annotate', 'spatial', 'cci', 'de', 'report']
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
    else if (params.stage == 'report')     rebuild_report()
    else                                   discover()
}

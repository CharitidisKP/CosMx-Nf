include { APPLY_LABELS } from '../modules/apply_labels'
include { SPATIAL_NET }  from '../modules/spatial_net'
include { BANKSY_PREP; BANKSY_SAMPLE; BANKSY } from '../modules/banksy'
include { BCELL }        from '../modules/bcell'
include { IMMUNE_DENSITY } from '../modules/immune_density'
include { COMPOSITION }  from '../modules/composition'
include { SMIDE }        from '../modules/smide'
include { SMIDE_META }   from '../modules/smide_meta'
include { SMIDE_OVERLAP } from '../modules/smide_overlap'
include { PSEUDOBULK_DE } from '../modules/pseudobulk_de'
include { BUNDLE }       from '../modules/bundle'
include { isIncluded }   from './samplesheet'

// Attaches cell_type from the run's label state, counts it and runs composition
// state: the basis, cluster labels, per-cell edits and digest from label_state
workflow ANNOTATION {
    take:
        merged_giotto
        state

    main:
        def st = [basis: state.basis, digest: state.digest,
                  labels_from: state.labels.toString(), overrides_from: state.overrides.toString()]
        ch = Channel.value( tuple('merged', file(merged_giotto), file(state.labels), file(state.overrides), st) )
        APPLY_LABELS(ch)
        COMPOSITION(APPLY_LABELS.out.obj, file(params.input))
        // The app bundle from the labelled object, its cells table passed on to REPORT
        def ch_bundle = channel.empty()
        if (params.bundle.toString().toBoolean()) {
            BUNDLE(APPLY_LABELS.out.obj)
            ch_bundle = BUNDLE.out.cells
        }

    emit:
        annotated = APPLY_LABELS.out.obj
        counts    = APPLY_LABELS.out.tables.flatten().filter { f -> f.name == 'cell_type_counts.csv' }
        tables    = APPLY_LABELS.out.tables.mix(COMPOSITION.out.tables, ch_bundle)
}

// The banksy_settings row, with the pipeline's lambdas and resolution
def banksy_settings() {
    def rows = file(params.banksy_variants).splitCsv(header: true)
    def row = rows.find { r -> r.variant == params.banksy_settings.toString() }
    if (!row)
        error "banksy_settings: no row '${params.banksy_settings}' in ${params.banksy_variants}. " +
              "Rows: ${rows.collect { r -> r.variant }.join(', ')}"
    row + [lambdas: params.banksy_lambda_ladder.toString(), resolutions: params.banksy_sample_res.toString()]
}

// What a spatial object's niches depend on besides the labels: a digest of the BANKSY settings
def banksy_key(v) {
    "${v.sort()}|${params.banksy_hvgs}|${params.banksy_bg_mult}|${params.banksy_match_cut}".toString().md5()
}

// Runs the Delaunay network, BANKSY, B-lineage subclustering and immune densities
// digest: the label state, saved beside the spatial object with the BANKSY digest
workflow SPATIAL {
    take:
        labelled
        digest

    main:
        SPATIAL_NET(labelled)
        def v = banksy_settings()
        BANKSY_PREP( SPATIAL_NET.out.obj.map { _id, obj -> obj }, '' )
        def ch_raw = BANKSY_PREP.out.raw.flatten()
            .map { f -> tuple(f.name.replaceFirst(/_raw\.rds$/, ''), f) }
        def ch_meta = BANKSY_PREP.out.meta.flatten()
            .map { f -> tuple(f.name.replaceFirst(/_meta\.csv\.gz$/, ''), f) }
        BANKSY_SAMPLE( ch_raw.join(ch_meta).map { s, r, m -> tuple(s, r, m, v) }, BANKSY_PREP.out.hvg.first() )
        BANKSY( SPATIAL_NET.out.obj, BANKSY_SAMPLE.out.labels.collect(), BANKSY_SAMPLE.out.genes.collect(),
                v, digest.map { d -> "${d}:${banksy_key(v)}".toString() } )
        BCELL(SPATIAL_NET.out.obj)
        IMMUNE_DENSITY(SPATIAL_NET.out.obj, file(params.input))

    emit:
        spatial = BANKSY.out.obj
        tables  = SPATIAL_NET.out.tables.mix(BANKSY.out.tables, BCELL.out.tables, IMMUNE_DENSITY.out.tables)
}

// Runs smiDE per sample and cell type, its meta-analysis, and pseudobulk DE
// counts: annotate's cell_type_counts.csv, read when smide_celltypes is empty
workflow DE {
    take:
        spatial
        counts

    main:
        ch_samples = Channel.fromPath(params.input).splitCsv(header: true)
                            .filter { isIncluded(it) }
                            .map { it.sample_id }.unique()
        def untyped = params.untyped_celltypes.toString().split(',').collect { it.trim() }
        ch_celltypes = params.smide_celltypes
            ? Channel.fromList(params.smide_celltypes.split(',').collect { it.trim() }.findAll { it })
            : counts.splitCsv(header: true)
                    .map { it.cell_type }.filter { it && !(it in untyped) }.unique()

        SMIDE( spatial.combine(ch_samples).combine(ch_celltypes) )
        // Overlap ratio per niche, beside the fits, so changing the flag never refits smiDE
        SMIDE_OVERLAP( spatial.combine(ch_samples), ch_celltypes.collect().map { it.join(',') } )
        SMIDE_META( SMIDE.out.de.collect(), SMIDE_OVERLAP.out.ovl.collect() )
        PSEUDOBULK_DE( spatial, file(params.input) )

    emit:
        tables = SMIDE_META.out.meta.mix(SMIDE_META.out.summary, PSEUDOBULK_DE.out.tables.flatten())
}

include { APPLY_LABELS } from '../modules/apply_labels'
include { SPATIAL_NET }  from '../modules/spatial_net'
include { BANKSY }       from '../modules/banksy'
include { BCELL }        from '../modules/bcell'
include { COMPOSITION }  from '../modules/composition'
include { SMIDE }        from '../modules/smide'
include { SMIDE_META }   from '../modules/smide_meta'
include { isIncluded }   from './samplesheet'

// Attach cell_type from the curated labels CSV, and count it. Nothing spatial.
workflow ANNOTATION {
    take:
        merged_giotto
        labels

    main:
        ch = Channel.value( tuple('merged', file(merged_giotto), file(labels)) )
        APPLY_LABELS(ch)
        COMPOSITION(APPLY_LABELS.out.obj, file(params.input))

    emit:
        annotated = APPLY_LABELS.out.obj
        tables    = APPLY_LABELS.out.tables.mix(COMPOSITION.out.tables)
}

// Delaunay network with proximity enrichment, then BANKSY domains.
// B-lineage subclustering needs neither, so it runs beside BANKSY.
workflow SPATIAL {
    take:
        labelled_giotto

    main:
        ch = Channel.value( tuple('merged', file(labelled_giotto)) )
        SPATIAL_NET(ch)
        BANKSY(SPATIAL_NET.out.obj)
        BCELL(SPATIAL_NET.out.obj)

    emit:
        spatial = BANKSY.out.obj
        tables  = SPATIAL_NET.out.tables.mix(BANKSY.out.tables, BCELL.out.tables)
}

// smiDE is per-sample by construction: the spatial model is per-tissue and the
// merged object is shift-joined. Fan out over (sample x cell type), gather by meta-analysis.
workflow DE {
    take:
        spatial_giotto
        labels

    main:
        ch = Channel.value( tuple('merged', file(spatial_giotto)) )
        ch_samples = Channel.fromPath(params.input).splitCsv(header: true)
                            .filter { isIncluded(it) }
                            .map { it.sample_id }.unique()
        ch_celltypes = params.smide_celltypes
            ? Channel.fromList(params.smide_celltypes.split(',').collect { it.trim() })
            : Channel.fromPath(labels).splitCsv(header: true)
                     .map { it.label }.filter { it && it.trim() }.unique()

        SMIDE( ch.combine(ch_samples).combine(ch_celltypes) )
        SMIDE_META( SMIDE.out.de.collect() )

    emit:
        tables = SMIDE_META.out.meta.mix(SMIDE_META.out.summary)
}

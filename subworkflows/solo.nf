include { INTEGRATE   as SOLO_INTEGRATE   } from '../modules/integrate'
include { CLUSTER     as SOLO_CLUSTER     } from '../modules/cluster'
include { HIERATYPE   as SOLO_HIERATYPE   } from '../modules/hieratype'
include { INSITUTYPE  as SOLO_INSITUTYPE  } from '../modules/insitutype'
include { DIAGNOSTICS as SOLO_DIAGNOSTICS } from '../modules/diagnostics'

// Discovery for one sample at a time, without MERGE or Harmony
workflow SOLO {
    take:
        qc_obj

    main:
        SOLO_INTEGRATE( qc_obj.map { meta, obj -> tuple(meta.id, obj) } )
        SOLO_CLUSTER( SOLO_INTEGRATE.out.obj )
        SOLO_HIERATYPE( SOLO_CLUSTER.out.obj )
        SOLO_INSITUTYPE( SOLO_HIERATYPE.out.obj )

        diag_cols = params.chosen_resolutions.toString().split(',').collect { "leiden_clus_res${it.trim()}" } + ['ht_call', 'insitutype_unsup']
        SOLO_DIAGNOSTICS( SOLO_INSITUTYPE.out.obj.combine(Channel.fromList(diag_cols)) )

    emit:
        obj = SOLO_INSITUTYPE.out.obj
}

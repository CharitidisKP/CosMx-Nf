include { CCI_LIANA }     from '../modules/cci_liana'
include { CCI_SPARKX }    from '../modules/cci_sparkx'
include { CCI_NNSVG }     from '../modules/cci_nnsvg'
include { CCI_INSITUCOR } from '../modules/cci_insitucor'
include { CCI_NICHENET }  from '../modules/cci_nichenet'
include { CCI_MISTY; CCI_MISTY_COLLECT } from '../modules/cci_misty'
include { isIncluded }    from './samplesheet'

workflow CCI {
    take:
        obj

    main:
        CCI_LIANA(obj)
        CCI_SPARKX(obj)
        CCI_INSITUCOR(obj)
        ch_samples = Channel.fromPath(params.input).splitCsv(header: true)
                            .filter { isIncluded(it) }
                            .map { it.sample_id }.unique()
        CCI_MISTY( obj.combine(ch_samples) )
        CCI_MISTY_COLLECT( CCI_MISTY.out.results.mix(CCI_MISTY.out.keys).collect() )
        ch_tables = CCI_LIANA.out.tables.mix(
            CCI_SPARKX.out.svgs, CCI_INSITUCOR.out.tables, CCI_MISTY_COLLECT.out.tables
        )
        // toBoolean, so --run_nnsvg false from the command line also skips
        if (params.run_nnsvg.toString().toBoolean()) {
            CCI_NNSVG(obj, CCI_SPARKX.out.svgs)
            ch_tables = ch_tables.mix(CCI_NNSVG.out.tables)
        }
        if (params.nichenet_model_dir) {
            CCI_NICHENET(obj, file(params.input))
            ch_tables = ch_tables.mix(CCI_NICHENET.out.tables)
        }

    emit:
        tables = ch_tables
}

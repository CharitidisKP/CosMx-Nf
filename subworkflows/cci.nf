include { CCI_LIANA }     from '../modules/cci_liana'
include { CCI_SPARKX }    from '../modules/cci_sparkx'
include { CCI_NNSVG }     from '../modules/cci_nnsvg'
include { CCI_INSITUCOR } from '../modules/cci_insitucor'
include { CCI_NICHENET }  from '../modules/cci_nichenet'
include { CCI_MISTY }     from '../modules/cci_misty'

workflow CCI {
    take:
        obj

    main:
        CCI_LIANA(obj)
        CCI_SPARKX(obj)
        CCI_NNSVG(obj, CCI_SPARKX.out.svgs)
        CCI_INSITUCOR(obj)
        CCI_MISTY(obj)
        ch_tables = CCI_LIANA.out.tables.mix(
            CCI_SPARKX.out.svgs, CCI_NNSVG.out.tables, CCI_INSITUCOR.out.tables, CCI_MISTY.out.tables
        )
        if (params.nichenet_model_dir) {
            CCI_NICHENET(obj, file(params.input))
            ch_tables = ch_tables.mix(CCI_NICHENET.out.tables)
        }

    emit:
        tables = ch_tables
}

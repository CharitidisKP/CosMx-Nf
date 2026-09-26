include { LOAD              } from '../modules/load'
include { QC                } from '../modules/qc'
include { MERGE             } from '../modules/merge'
include { INTEGRATE         } from '../modules/integrate'
include { CLUSTER           } from '../modules/cluster'
include { HIERATYPE         } from '../modules/hieratype'
include { INSITUTYPE        } from '../modules/insitutype'
include { SUPERVISED ; ATTACH_SUPERVISED } from '../modules/insitutype_supervised'
include { DIAGNOSTICS       } from '../modules/diagnostics'
include { SOLO              } from './solo'
include { isIncluded        } from './samplesheet'

workflow DISCOVERY {
    take:
    samplesheet
    solo_ids

    main:
    ch_samples = Channel
        .fromPath(samplesheet)
        .splitCsv(header: true)
        .filter { row -> isIncluded(row) }
        .filter { row -> !params.solo_only || solo_ids.contains(row.sample_id) }
        .map { row ->
            tuple(
                [id: row.sample_id, prefix: row.file_prefix, patient: row.patient_id, treatment: row.treatment, timepoint: row.timepoint, slide: row.slide_id, batch: row.batch, subset_fovs: row.subset_fovs],
                file(row.data_dir),
            )
        }

    LOAD(ch_samples)
    QC(LOAD.out.obj)

    ch_merged = Channel.empty()
    ch_tables = LOAD.out.tables.mix(QC.out.tables)

    if (!params.solo_only) {
        ch_merge = QC.out.obj
            .toSortedList { a, b -> a[0].id <=> b[0].id }
            .map { rows -> tuple(rows.collect { it[0].id }, rows.collect { it[1] }) }

        MERGE(ch_merge)
        INTEGRATE(MERGE.out.obj)
        CLUSTER(INTEGRATE.out.obj)
        HIERATYPE(CLUSTER.out.obj)
        INSITUTYPE(HIERATYPE.out.obj)

        // Profiles are named in the config and staged per machine under profile_dir.
        def prof_names = (params.supervised_profiles ?: "")
            .toString()
            .split(',')
            .collect { n -> n.trim() }
            .findAll { n -> n }
        def profs = prof_names.collect { n ->
            [name: n, reference: ['rds', 'csv'].collect { ext -> file("${params.profile_dir}/${n}.${ext}") }.find { f -> f.exists() }]
        }
        def missing = profs.findAll { p -> !p.reference }
        if (missing) {
            error(
                "supervised_profiles: no .rds or .csv in ${params.profile_dir} for:\n" + missing.collect { p -> "  ${p.name}" }.join('\n') + "\nStage the profile matrices there, or set supervised_profiles = '' to skip."
            )
        }

        ch_tables = ch_tables.mix(
            MERGE.out.tables, INTEGRATE.out.tables, CLUSTER.out.tables,
            HIERATYPE.out.tables, INSITUTYPE.out.tables
        )

        if (profs) {
            SUPERVISED(INSITUTYPE.out.obj.combine(Channel.fromList(profs)))
            ch_tables = ch_tables.mix(SUPERVISED.out.tables)
            ch_sup = SUPERVISED.out.cells
                .map { pid, csv -> tuple(pid, csv) }
                .groupTuple(size: profs.size())
            ATTACH_SUPERVISED(INSITUTYPE.out.obj.join(ch_sup))
            ch_merged = ATTACH_SUPERVISED.out.obj
        }
        else {
            ch_merged = INSITUTYPE.out.obj
        }

        def sup_cols = profs
            .collect { p ->
                params.insitutype_refine
                    ? ["insitutype_sup_${p.name}", "insitutype_sup_${p.name}_refined"]
                    : ["insitutype_sup_${p.name}"]
            }
            .flatten()
        diag_cols = params.chosen_resolutions.toString().split(',').collect { "leiden_clus_res${it.trim()}" } + ['ht_call', 'insitutype_unsup'] + sup_cols
        DIAGNOSTICS(ch_merged.combine(Channel.fromList(diag_cols)))
        ch_tables = ch_tables.mix(DIAGNOSTICS.out.tables)
    }

    // Opt-in: take named samples through integrate -> cluster -> insitutype -> diagnostics on their own.
    if (solo_ids) {
        SOLO(QC.out.obj.filter { meta, obj -> solo_ids.contains(meta.id) })
    }

    emit:
    merged = ch_merged
    tables = ch_tables
}

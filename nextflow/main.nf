#!/usr/bin/env nextflow

// Disease/BMI run, adapted from Koslicki_R1/nextflow. Differences from Koslicki_R1:
//   * input is a todo TSV (scripts/make_todo.py) that says, per accession, which stages still
//     have to run -- anything already in the shared all-res tree is not recomputed;
//   * accessions that only lack YACHT go straight to YACHT from their existing all-res sketch;
//   * results publish into all-res (params.outdir), so its completed-accession index stays current.
//
// Usage: nextflow run nextflow --todo shard_000.tsv [-resume]

nextflow.enable.dsl = 2

include { FETCH_LOGAN }      from './modules/fetch.nf'
include { WRITE_LEDGER_ROW } from './modules/ledger.nf'
include { ANALYZE }          from './subworkflows/analyze.nf'

workflow {

    todo = Channel.fromPath(params.todo, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { r ->
            def need = [
                unitigs: [sketch: r.sketch_unitigs == '1', ko: r.ko_unitigs == '1', yacht: r.yacht_unitigs == '1'],
                contigs: [sketch: r.sketch_contigs == '1', ko: r.ko_contigs == '1', yacht: false],
            ]
            [r.accession, r.fetch_unitigs == '1', r.fetch_contigs == '1', need, r.existing_unitigs_sig]
        }

    // needs: accession -> per-seq_type stage flags, re-attached to fetch outputs below
    needs = todo.map { acc, fu, fc, need, sig -> tuple(acc, need) }

    fetched = FETCH_LOGAN(todo.filter { it[1] || it[2] }.map { acc, fu, fc, need, sig -> tuple(acc, fu, fc) })
    fetch_ch = fetched.unitigs.mix(fetched.contigs)   // tuple(accession, seq_type, zst_path, status)

    fetch_ch.branch {
        ok:    it[3] == 'OK'
        empty: it[3] == 'SKIPPED_EMPTY'
        // NOT_REQUESTED: this seq_type's outputs already exist in all-res
    }.set { gated }

    WRITE_LEDGER_ROW(
        gated.empty.map { row -> tuple(row[0], row[1], 'fetch', row[3]) }
    )

    to_analyze = gated.ok
        .map { acc, seq_type, zst, status -> tuple(acc, seq_type, zst) }
        .combine(needs, by: 0)
        .map { acc, seq_type, zst, need -> tuple(acc, seq_type, zst, need[seq_type]) }

    existing_sigs = todo
        .filter { acc, fu, fc, need, sig -> sig }
        .map { acc, fu, fc, need, sig -> tuple(acc, 'unitigs', file(sig, checkIfExists: true)) }

    ANALYZE(to_analyze, existing_sigs)
}

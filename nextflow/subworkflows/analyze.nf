// Analysis fan-out, from Koslicki_R1/nextflow/subworkflows/analyze.nf, with per-stage gating:
// each decompressed FASTA only goes to the stages whose output is still missing from all-res
// (need.sketch / need.ko / need.yacht, from scripts/make_todo.py). YACHT runs on unitigs only,
// on either a freshly made sketch or one that already exists in all-res.

include { DECOMPRESS_FASTA } from '../modules/decompress.nf'
include { SOURMASH_SKETCH }  from '../modules/sketch.nf'
include { FUNPROFILER }      from '../modules/funprofiler.nf'
include { YACHT }            from '../modules/yacht.nf'
include { CLEANUP_FASTA }    from '../modules/cleanup.nf'

workflow ANALYZE {
    take:
    zst_fasta_ch    // tuple(accession, seq_type, zst_fasta_path, need)  need = [sketch:, ko:, yacht:]
    existing_sigs   // tuple(accession, 'unitigs', sig_zip) -- YACHT-only accessions

    main:
    // keep the stage flags keyed by accession.seq_type; they can't travel through the processes
    need_by_key = zst_fasta_ch.map { acc, seq_type, zst, need -> tuple("${acc}.${seq_type}", need) }

    decompressed = DECOMPRESS_FASTA(zst_fasta_ch.map { acc, seq_type, zst, need -> tuple(acc, seq_type, zst) })
    fasta_need = decompressed.fasta
        .map { acc, seq_type, fasta -> tuple("${acc}.${seq_type}", acc, seq_type, fasta) }
        .join(need_by_key)          // tuple(key, acc, seq_type, fasta, need)

    sketch_out = SOURMASH_SKETCH(fasta_need.filter { it[4].sketch }.map { k, acc, st, fa, need -> tuple(acc, st, fa) })
    ko_out     = FUNPROFILER(fasta_need.filter { it[4].ko }.map { k, acc, st, fa, need -> tuple(acc, st, fa) })

    // a new sketch goes to YACHT only if that accession still lacks a YACHT result
    new_sigs_for_yacht = sketch_out.sig
        .filter { acc, seq_type, sig -> seq_type == 'unitigs' }
        .map { acc, seq_type, sig -> tuple("${acc}.${seq_type}", acc, seq_type, sig) }
        .join(need_by_key)
        .filter { it[4].yacht }
        .map { k, acc, st, sig, need -> tuple(acc, st, sig) }
    yacht_out = YACHT(new_sigs_for_yacht.mix(existing_sigs))

    // Delete the (often multi-GB) FASTA as soon as every stage that needs it is done; a stage that
    // isn't needed counts as done immediately.
    sketch_done = sketch_out.sig.map { acc, st, sig -> tuple("${acc}.${st}", true) }
        .mix(fasta_need.filter { !it[4].sketch }.map { tuple(it[0], true) })
    ko_done = ko_out.ko_csv.map { acc, st, csv -> tuple("${acc}.${st}", true) }
        .mix(fasta_need.filter { !it[4].ko }.map { tuple(it[0], true) })

    CLEANUP_FASTA(
        fasta_need.map { k, acc, st, fa, need -> tuple(k, acc, st, fa) }
            .join(sketch_done).join(ko_done)
            .map { k, acc, st, fa, _s, _k -> tuple(acc, st, fa) }
    )

    emit:
    sig    = sketch_out.sig
    ko     = ko_out.ko_csv
    xlsx   = yacht_out.xlsx
    ledger = sketch_out.ledger.mix(ko_out.ledger, yacht_out.ledger)
}

// YACHT ANI-based presence/absence, matching the yacht invocation in
// logan-subworkflow.sh. `params.yacht_ani_thresh` is only used for the output filename
// (mirroring the reference script) -- it must match the `ani_thresh` field already
// baked into the pretrained manifest referenced by `params.yacht_config`.

process YACHT {
    tag "${accession}.${seq_type}"
    // path must be a closure -- see the NOTE in modules/sketch.nf.
    publishDir(path: { "${params.outdir}/yacht_results/${seq_type}" }, mode: 'copy', pattern: "*.xlsx")
    publishDir "${params.outdir}/ledger", mode: 'copy', pattern: "*.ledger.csv"

    // YACHT's RAM grows with how many reference genomes overlap the sample, which tracks the
    // sketch size. The local executor never enforces this -- it only decides how many tasks
    // start (see the executor{} cap in nextflow.config), so it should match real usage.
    // History: 8/16/24/32/40/48 GB on n=1000 kept the 2400 GB cap full while the server's actual
    // RAM use was only ~52% (avg) / ~60% (peak) of that, with no swap -- scaled to 62% here.
    // Sketches over 8 MB were rare in those runs (<10% of tasks) so a 30 GB top tier averaged
    // out; for batches of large metagenomes it would under-reserve. Nextflow's peak_rss for
    // 8-20 MB sketches was ~80 GB median, ~1.8x real usage -> ~45 GB. Above 20 MB is unmeasured:
    // no benchmark included an accession over the old 2 GB sample filter. Sketch size scales
    // ~linearly with the compressed unitigs file (~15 MB at 2 GB, ~67 MB at 10 GB, ~230 MB for the
    // largest 37.6 GB accession), so unfiltered samples reach well past 20 MB; the 80 GB tier is a
    // guess for them -- check the RAM panel of plot_benchmark.py for swap after such runs.
    memory {
        def mb = sig_zip.size() / 1_000_000
        mb < 0.1 ? 5.GB : mb < 0.5 ? 10.GB : mb < 1 ? 15.GB : mb < 2 ? 20.GB : mb < 4 ? 25.GB :
            mb < 8 ? 30.GB : mb < 16 ? 45.GB : mb < 50 ? 60.GB : 80.GB
    }

    input:
    tuple val(accession), val(seq_type), path(sig_zip)

    output:
    tuple val(accession), val(seq_type), path("result_${params.yacht_ani_thresh}_ANI_${accession}_${seq_type}.xlsx"), emit: xlsx
    path("${accession}.${seq_type}.yacht.ledger.csv"), emit: ledger

    script:
    """
    set -euo pipefail

    out="result_${params.yacht_ani_thresh}_ANI_${accession}_${seq_type}.xlsx"

    # A tiny assembly can yield a sketch with zero hashes at scaled=1000; YACHT then errors with
    # "Empty sketch in signature". That is a valid nothing-to-profile outcome, not a failure.
    n_hashes=\$(python -c "import csv,sys,zipfile; z=zipfile.ZipFile(sys.argv[1]); r=list(csv.DictReader(l for l in z.read('SOURMASH-MANIFEST.csv').decode().splitlines() if not l.startswith('#'))); print(sum(int(x['n_hashes']) for x in r if x['ksize']=='31'))" "${sig_zip}")
    if [ "\$n_hashes" -eq 0 ]; then
        : > "\$out"
        write_ledger_row.sh "${accession}.${seq_type}.yacht.ledger.csv" "${accession}" "${seq_type}" yacht SKIPPED_EMPTY_SKETCH
        exit 0
    fi

    set +e
    yacht run --json "${params.yacht_config}" \\
        --sample_file "${sig_zip}" --significance ${params.yacht_significance} \\
        --num_threads ${task.cpus} \\
        --min_coverage_list ${params.yacht_min_coverage_list} \\
        --out "\$out"
    yacht_rc=\$?
    set -e

    # If YACHT returned a nonzero exit code, the script exits with the same code. Nextflow then marks the YACHT task as failed.
    if [ \$yacht_rc -ne 0 ]; then
        exit \$yacht_rc
    fi

    # check if yacht found any hits
    if [ -s "\$out" ]; then
        write_ledger_row.sh "${accession}.${seq_type}.yacht.ledger.csv" "${accession}" "${seq_type}" yacht DONE
    else
        : > "\$out"
        write_ledger_row.sh "${accession}.${seq_type}.yacht.ledger.csv" "${accession}" "${seq_type}" yacht SKIPPED_NO_HITS
    fi
    """
}

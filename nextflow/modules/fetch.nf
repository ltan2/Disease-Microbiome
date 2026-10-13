// Fetch + empty gate, from Koslicki_R1/nextflow/modules/fetch.nf. Only the sequence types that
// still need a sketch or KO profile are downloaded; the other one is emitted as NOT_REQUESTED
// (an empty placeholder file) so both outputs always exist.

process FETCH_LOGAN {
    tag "$accession"

    // Transient S3/network failures are retried; once retries are exhausted the accession is
    // skipped (params.fetch_exhausted_strategy) rather than aborting a multi-week run -- it will
    // show as incomplete in all-res/index and is picked up by the next make_todo.py.
    errorStrategy { task.attempt <= task.maxRetries ? 'retry' : params.fetch_exhausted_strategy }
    maxRetries 3

    input:
    tuple val(accession), val(fetch_unitigs), val(fetch_contigs)

    output:
    tuple val(accession), val('unitigs'), path("${accession}.unitigs.fa.zst"), env('U_STATUS'), emit: unitigs
    tuple val(accession), val('contigs'), path("${accession}.contigs.fa.zst"), env('C_STATUS'), emit: contigs

    script:
    """
    set -euo pipefail

    # Backpressure (replaces sharding): wait while too many fetched sequences are still waiting
    # for sketch/funprofiler. One token per OK sequence lives in params.inflight_dir from here until
    # CLEANUP_FASTA deletes its FASTA; tokens older than a day (a failed downstream task) are ignored.
    mkdir -p "${params.inflight_dir}"
    while [ "\$(find "${params.inflight_dir}" -type f -mmin -1440 | wc -l)" -ge ${params.max_inflight} ]; do
        sleep 20
    done

    # Logan doesn't always publish both files for a given accession; a missing key is a
    # permanent 404, so treat it as an empty download (caught by the size gate below).
    fetch_or_empty() {
        local uri="\$1" out="\$2"
        if aws s3 cp "\$uri" "\$out" --no-sign-request --only-show-errors 2>fetch.err; then
            return 0
        elif grep -q "does not exist\\|(404)" fetch.err; then
            : > "\$out"
            return 0
        else
            cat fetch.err >&2
            return 1
        fi
    }

    # size gate: zstd hangs on near-empty Logan files (reference-logan-analysis.sh)
    gate() {
        if [ "\$(stat -c%s "\$1")" -lt ${params.min_compressed_bytes} ]; then echo SKIPPED_EMPTY; else echo OK; fi
    }

    if [ "${fetch_unitigs}" = "true" ]; then
        fetch_or_empty "s3://logan-pub/u/${accession}/${accession}.unitigs.fa.zst" "${accession}.unitigs.fa.zst"
        U_STATUS=\$(gate "${accession}.unitigs.fa.zst")
        if [ "\$U_STATUS" = OK ]; then touch "${params.inflight_dir}/${accession}.unitigs"; fi
    else
        : > "${accession}.unitigs.fa.zst"
        U_STATUS=NOT_REQUESTED
    fi

    if [ "${fetch_contigs}" = "true" ]; then
        fetch_or_empty "s3://logan-pub/c/${accession}/${accession}.contigs.fa.zst" "${accession}.contigs.fa.zst"
        C_STATUS=\$(gate "${accession}.contigs.fa.zst")
        if [ "\$C_STATUS" = OK ]; then touch "${params.inflight_dir}/${accession}.contigs"; fi
    else
        : > "${accession}.contigs.fa.zst"
        C_STATUS=NOT_REQUESTED
    fi
    """
}

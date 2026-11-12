// Stages "merge" and "combine": consensus calls, orfipy novel ORFs, hygiene.

process MERGE_CHUNK {
    tag "${chunk_id}"
    label 'py_chunk'

    input:
    tuple val(chunk_id), path(genomes), path(stats), path(bakta_dirs), path(smorf_dirs), path(gmsc_dir), path(extra_dirs)
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    tuple val(chunk_id), path("${chunk_id}.merged", type: 'dir'), path("${chunk_id}.merged_summary.tsv"), emit: merged
    tuple val(chunk_id), path("${chunk_id}.merged_orf_table.tsv"), path("${chunk_id}.merged_summary.tsv"), emit: table
    tuple val(chunk_id), path("${chunk_id}.merge_skipped.tsv"), path("${chunk_id}.merge_combos.tsv"), emit: info

    script:
    // Extra callers merge after the built-in three, in params order (ties keep the earlier call).
    def extra = params.extra_callers
        ? "--extra-tools ${params.extra_callers.collect { c -> c.name }.join(' ')} --extra-dirs ${extra_dirs}" : ''
    """
    ${params.python_bin} ${projectDir}/bin/merge_chunk.py \\
        --genome-stats ${stats} \\
        --genomes-dir ${genomes} \\
        --bakta-dirs ${bakta_dirs} \\
        --smorf-dirs ${smorf_dirs} \\
        --gmsc-dir ${gmsc_dir} \\
        --max-aa ${params.max_aa} \\
        --prefix ${chunk_id} ${extra}
    """

    stub:
    """
    mkdir ${chunk_id}.merged
    touch ${chunk_id}.merged_summary.tsv ${chunk_id}.merged_orf_table.tsv ${chunk_id}.merge_skipped.tsv ${chunk_id}.merge_combos.tsv
    """
}

process ORFIPY_NOVEL {
    tag "${chunk_id}"
    label 'orfipy'

    input:
    tuple val(chunk_id), path(genomes), path(merged_dir), path(merged_summary)
    path code, stageAs: 'code/*'

    output:
    tuple val(chunk_id), path("${chunk_id}.orfipy", type: 'dir'), path("${chunk_id}.orfipy_stats.tsv")

    script:
    // orfipy_worker.sh:30-34,123-129: exhaustive ORFs 15-300 nt with the stop codon included,
    // then filter_novel_orfs.py drops any ORF touching an annotated feature (either strand).
    // Each genome is copied in first: orfipy writes an index next to its input FASTA.
    """
    mkdir -p ${chunk_id}.orfipy
    printf 'genome_id\\tn_raw\\tn_novel\\tn_gene_features\\tstatus\\n' > ${chunk_id}.orfipy_stats.tsv
    while read -r gid; do
        w="w_\$gid"; rm -rf "\$w"; mkdir -p "\$w/out"
        cp "${genomes}/\$gid.fna" "\$w/\$gid.fna"
        rc=0
        timeout 20m ${params.orfipy_bin} \\
            --min ${params.orf_min} --max ${params.orf_max} \\
            --start ${params.orf_starts} --stop ${params.orf_stops} \\
            --strand b --include-stop \\
            --procs ${task.cpus} --chunk-size 100 \\
            --outdir "\$w/out" --bed "\$gid.bed" \\
            "\$w/\$gid.fna" > "\$w/orfipy.log" 2>&1 || rc=\$?
        if [ "\$rc" -ne 0 ] || [ ! -f "\$w/out/\$gid.bed" ]; then
            echo "\$gid orfipy failed rc=\$rc" >&2; tail -5 "\$w/orfipy.log" >&2 || true
            printf '%s\\t\\t\\t\\torfipy_rc_%s\\n' "\$gid" "\$rc" >> ${chunk_id}.orfipy_stats.tsv
            rm -rf "\$w"; continue
        fi
        ${params.python_bin} ${projectDir}/bin/filter_novel_orfs.py \\
            --bed "\$w/out/\$gid.bed" \\
            --merged-gff "${merged_dir}/\$gid.gff" \\
            --genome-id "\$gid" \\
            --out "${chunk_id}.orfipy/\$gid.gff" \\
            --stats "\$w/stats.tsv"
        printf '%s\\tok\\n' "\$(cat "\$w/stats.tsv")" >> ${chunk_id}.orfipy_stats.tsv
        rm -rf "\$w"
    done < <(tail -n +2 ${merged_summary} | cut -f1)
    """

    stub:
    """
    mkdir ${chunk_id}.orfipy
    touch ${chunk_id}.orfipy_stats.tsv
    """
}

process COMBINE_CHUNK {
    tag "${chunk_id}"
    label 'py_chunk'

    input:
    tuple val(chunk_id), path(genomes), path(merged_table), path(merged_summary), path(orfipy_dir), path(orfipy_stats)
    path hygiene
    path code, stageAs: 'code/*'

    output:
    tuple val(chunk_id), path("${chunk_id}.smorfs.faa"), path("${chunk_id}.smorfs_metadata.tsv"), emit: smorfs
    tuple val(chunk_id), path("${chunk_id}.dropped.tsv"), path("${chunk_id}.combine_map.tsv"), path("${chunk_id}.combine_stats.tsv"), emit: info

    script:
    """
    ${params.python_bin} ${projectDir}/bin/combine_chunk.py \\
        --merged-table ${merged_table} \\
        --merged-summary ${merged_summary} \\
        --orfipy-dir ${orfipy_dir} \\
        --genomes-dir ${genomes} \\
        --hygiene ${hygiene} \\
        --max-aa ${params.max_aa} \\
        --prefix ${chunk_id}
    """

    stub:
    """
    touch ${chunk_id}.smorfs.faa ${chunk_id}.smorfs_metadata.tsv ${chunk_id}.dropped.tsv ${chunk_id}.combine_map.tsv ${chunk_id}.combine_stats.tsv
    """
}

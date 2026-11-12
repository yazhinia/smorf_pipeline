// Stage "flags": AntiFam and Popcorn on each distinct sequence, joined back to every smORF.

process EXTRACT_SEQS {
    tag "${chunk_id}"
    label 'py_chunk'

    input:
    tuple val(chunk_id), path(genomes), path(faa), path(metadata)
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    tuple val(chunk_id), path("${chunk_id}.aa_unique.faa"), emit: aa
    tuple val(chunk_id), path("${chunk_id}.nt_unique.ffn"), emit: nt
    tuple val(chunk_id), path("${chunk_id}.nt_map.tsv"), emit: nt_map

    script:
    """
    ${params.python_bin} ${projectDir}/bin/extract_seqs.py \\
        --faa ${faa} --metadata ${metadata} --genomes-dir ${genomes} --prefix ${chunk_id}
    """

    stub:
    """
    touch ${chunk_id}.aa_unique.faa ${chunk_id}.nt_unique.ffn ${chunk_id}.nt_map.tsv
    """
}

process ANTIFAM {
    tag "${chunk_id}"
    label 'hmmer'

    input:
    tuple val(chunk_id), path(aa)

    output:
    tuple val(chunk_id), path("${chunk_id}.antifam.tbl")

    script:
    // antifam_worker.sh:98-114. -Z pins E-values to the catalog size; --cut_ga thresholds on
    // bit score, so -Z never changes which sequences are hit.
    """
    source ${projectDir}/bin/pam_modules.sh ${params.hmmer_module}
    command -v hmmsearch > /dev/null || { echo "hmmsearch not on PATH after module load" >&2; exit 3; }
    if [ -s ${aa} ]; then
        hmmsearch --cut_ga --cpu ${task.cpus} -Z ${params.antifam_z} \\
            -o /dev/null --tblout ${chunk_id}.antifam.tbl ${params.antifam_hmm} ${aa}
    else
        printf '# no sequences\\n# [ok]\\n' > ${chunk_id}.antifam.tbl
    fi
    tail -1 ${chunk_id}.antifam.tbl | grep -q '^# \\[ok\\]' || { echo "hmmsearch output incomplete" >&2; exit 4; }
    """

    stub:
    """
    printf '# [ok]\\n' > ${chunk_id}.antifam.tbl
    """
}

process POPCORN {
    tag "${chunk_id}"
    label 'popcorn'

    input:
    tuple val(chunk_id), path(nt)

    output:
    tuple val(chunk_id), path("${chunk_id}.popcorn.tsv")

    script:
    // run_popcorn.sh:56-59: generic model, -m/-z passed explicitly so the working
    // directory does not matter.
    """
    if [ -s ${nt} ]; then
        ${params.popcorn_python} ${params.popcorn_dir}/Popcorn.py \\
            -f ${nt} -o ${chunk_id}.popcorn.tsv \\
            -m ${params.popcorn_dir}/DataFiles/model.all.pickle \\
            -z ${params.popcorn_dir}/DataFiles/genome.all.pickle
    else
        printf 'Sequence ID\\tCoding Probability\\tPrediction\\n' > ${chunk_id}.popcorn.tsv
    fi
    """

    stub:
    """
    printf 'Sequence ID\\tCoding Probability\\tPrediction\\n' > ${chunk_id}.popcorn.tsv
    """
}

process JOIN_FLAGS {
    tag "${chunk_id}"
    label 'py_chunk'

    input:
    tuple val(chunk_id), path(faa), path(metadata), path(nt_map), path(tbl), path(popcorn), path(contig_map)
    path code, stageAs: 'code/*'

    output:
    tuple val(chunk_id), path("${chunk_id}.smorfs.tsv"), emit: smorfs
    tuple val(chunk_id), path("${chunk_id}.flags_stats.tsv"), emit: stats

    script:
    """
    ${params.python_bin} ${projectDir}/bin/join_flags.py \\
        --metadata ${metadata} --faa ${faa} --nt-map ${nt_map} \\
        --antifam ${tbl} --popcorn ${popcorn} --contig-map ${contig_map} \\
        --prefix ${chunk_id}
    """

    stub:
    """
    touch ${chunk_id}.smorfs.tsv ${chunk_id}.flags_stats.tsv
    """
}

// Checkpoint summaries and the final run-level tables.

process STAGE_SUMMARY {
    tag "${stage}"
    label 'py_small'

    input:
    tuple val(stage), path(inputs)
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    path 'qc/*'

    script:
    """
    ${params.python_bin} ${projectDir}/bin/stage_summary.py --stage ${stage} --out-dir qc
    """

    stub:
    """
    mkdir qc
    printf '## Stage `${stage}`: OK (stub)\\n' > qc/${stage}_summary.md
    printf 'metric\\tvalue\\tlevel\\tnote\\n' > qc/${stage}_summary.tsv
    """
}

process COLLECT_RESULTS {
    label 'py_chunk'

    input:
    path inputs
    path code, stageAs: 'code/*'

    output:
    path '*.{tsv,faa,gff3}'

    script:
    """
    ${params.python_bin} ${projectDir}/bin/collect_results.py \\
        --concat-max-genomes ${params.concat_max_genomes}
    """

    stub:
    """
    touch per_genome_summary.tsv genome_map.tsv smorfs.tsv
    """
}

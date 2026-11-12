// Stage "prepare": decide the genomes, then normalize them chunk by chunk.

process PLAN_INPUTS {
    label 'py_small'

    input:
    val input_path
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    path 'plan.tsv', emit: plan

    script:
    def regex = params.split_regex ? "--split-regex '${params.split_regex}'" : ''
    """
    ${params.python_bin} ${projectDir}/bin/plan_inputs.py \\
        --input '${input_path}' \\
        --split-mode ${params.split_mode} ${regex} \\
        --id-prefix ${params.id_prefix} \\
        --chunk-size ${params.chunk_size} \\
        --out plan.tsv
    """

    stub:
    """
    printf 'genome_id\\tsample_name\\tsource_fasta\\tselector\\tchunk_id\\n' > plan.tsv
    printf 'gU000001\\ts1\\t${input_path}\\t*\\tchunk_000001\\n' >> plan.tsv
    printf 'gU000002\\ts2\\t${input_path}\\t*\\tchunk_000001\\n' >> plan.tsv
    """
}

process NORMALIZE_CHUNK {
    tag "${chunk_id}"
    label 'py_chunk'

    input:
    tuple val(chunk_id), path(plan)
    path code, stageAs: 'code/*'

    output:
    tuple val(chunk_id), path("${chunk_id}.genomes", type: 'dir'), emit: genomes
    tuple val(chunk_id), path("${chunk_id}.fna"), path("${chunk_id}.contig_map.tsv"), emit: chunk_fasta
    tuple val(chunk_id), path("${chunk_id}.genome_stats.tsv"), emit: stats
    tuple val(chunk_id), path("${chunk_id}.contig_map.tsv"), emit: contig_map

    script:
    def regex = params.split_regex ? "--split-regex '${params.split_regex}'" : ''
    """
    ${params.python_bin} ${projectDir}/bin/normalize_chunk.py \\
        --plan ${plan} \\
        --chunk-id ${chunk_id} \\
        --split-mode ${params.split_mode} ${regex} \\
        --prefix ${chunk_id}
    """

    stub:
    """
    mkdir ${chunk_id}.genomes
    printf '>gU000001_c001\\nATG\\n' > ${chunk_id}.genomes/gU000001.fna
    printf '>gU000002_c001\\nATG\\n' > ${chunk_id}.genomes/gU000002.fna
    cat ${chunk_id}.genomes/*.fna > ${chunk_id}.fna
    touch ${chunk_id}.contig_map.tsv ${chunk_id}.genome_stats.tsv
    """
}

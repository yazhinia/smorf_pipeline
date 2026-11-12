include { PLAN_INPUTS; NORMALIZE_CHUNK } from '../../modules/local/prepare'
include { STAGE_SUMMARY } from '../../modules/local/report'
include { binFiles } from '../../modules/local/utils'

workflow PREPARE {
    take:
    input_path

    main:
    PLAN_INPUTS(input_path, binFiles(['plan_inputs.py']))

    // Chunk ids are 'chunk_000001'..; --max_chunks keeps the first N for a pilot run.
    // The plan is deterministic, so a later full run reuses these chunks' tasks.
    chunks = PLAN_INPUTS.out.plan
        .splitCsv(header: true, sep: '\t')
        .map { row -> row.chunk_id }
        .unique()
        .filter { cid -> params.max_chunks <= 0 || (cid.tokenize('_')[1] as Integer) <= params.max_chunks }

    NORMALIZE_CHUNK(chunks.combine(PLAN_INPUTS.out.plan), binFiles(['normalize_chunk.py']))

    STAGE_SUMMARY(
        NORMALIZE_CHUNK.out.stats.map { _cid, f -> f }.collect().map { fs -> tuple('prepare', fs) },
        binFiles(['stage_summary.py'])
    )

    emit:
    genomes     = NORMALIZE_CHUNK.out.genomes
    chunk_fasta = NORMALIZE_CHUNK.out.chunk_fasta
    stats       = NORMALIZE_CHUNK.out.stats
    contig_map  = NORMALIZE_CHUNK.out.contig_map
}

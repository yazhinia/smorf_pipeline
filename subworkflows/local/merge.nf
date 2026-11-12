include { MERGE_CHUNK } from '../../modules/local/merge_combine'
include { STAGE_SUMMARY } from '../../modules/local/report'
include { binFiles } from '../../modules/local/utils'

workflow MERGE {
    take:
    genomes    // tuple(chunk_id, genomes_dir)
    stats      // tuple(chunk_id, genome_stats)
    bakta      // tuple(chunk_id, [dirs])
    smorf      // tuple(chunk_id, [dirs])
    gmsc       // tuple(chunk_id, dir)
    extra      // tuple(chunk_id, [dirs]), [] without extra_callers

    main:
    MERGE_CHUNK(genomes.join(stats).join(bakta).join(smorf).join(gmsc).join(extra),
                binFiles(['merge_chunk.py', 'merge_gffs.py']))

    summary_files = MERGE_CHUNK.out.table.map { _cid, _t, s -> s }
        .mix(MERGE_CHUNK.out.info.flatMap { _cid, skipped, combos -> [skipped, combos] })
        .collect()
    STAGE_SUMMARY(summary_files.map { fs -> tuple('merge', fs) }, binFiles(['stage_summary.py']))

    emit:
    merged = MERGE_CHUNK.out.merged   // tuple(chunk_id, merged_dir, merged_summary)
    table  = MERGE_CHUNK.out.table    // tuple(chunk_id, merged_orf_table, merged_summary)
    info   = MERGE_CHUNK.out.info     // tuple(chunk_id, skipped, combos)
}

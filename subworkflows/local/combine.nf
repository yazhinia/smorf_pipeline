include { ORFIPY_NOVEL; COMBINE_CHUNK } from '../../modules/local/merge_combine'
include { STAGE_SUMMARY } from '../../modules/local/report'
include { binFiles } from '../../modules/local/utils'

workflow COMBINE {
    take:
    genomes    // tuple(chunk_id, genomes_dir)
    merged     // tuple(chunk_id, merged_dir, merged_summary)
    table      // tuple(chunk_id, merged_orf_table, merged_summary)

    main:
    ORFIPY_NOVEL(genomes.join(merged), binFiles(['filter_novel_orfs.py', 'merge_gffs.py']))
    COMBINE_CHUNK(
        genomes.join(table).join(ORFIPY_NOVEL.out),
        file("${projectDir}/assets/hygiene.yaml", checkIfExists: true),
        binFiles(['combine_chunk.py', 'hygiene.py', 'merge_chunk.py', 'merge_gffs.py'])
    )

    summary_files = COMBINE_CHUNK.out.info.flatMap { _cid, _dropped, cmap, cstats -> [cmap, cstats] }
        .collect()
    STAGE_SUMMARY(summary_files.map { fs -> tuple('combine', fs) }, binFiles(['stage_summary.py']))

    emit:
    smorfs = COMBINE_CHUNK.out.smorfs   // tuple(chunk_id, faa, metadata)
    info   = COMBINE_CHUNK.out.info     // tuple(chunk_id, dropped, combine_map, combine_stats)
}

include { EXTRACT_SEQS; ANTIFAM; POPCORN; JOIN_FLAGS } from '../../modules/local/flags'
include { STAGE_SUMMARY } from '../../modules/local/report'
include { binFiles } from '../../modules/local/utils'

workflow FLAGS {
    take:
    genomes      // tuple(chunk_id, genomes_dir)
    smorfs       // tuple(chunk_id, faa, metadata)
    contig_map   // tuple(chunk_id, contig_map)

    main:
    EXTRACT_SEQS(genomes.join(smorfs), binFiles(['extract_seqs.py', 'merge_chunk.py', 'merge_gffs.py']))
    ANTIFAM(EXTRACT_SEQS.out.aa)
    POPCORN(EXTRACT_SEQS.out.nt)
    JOIN_FLAGS(
        smorfs.join(EXTRACT_SEQS.out.nt_map).join(ANTIFAM.out).join(POPCORN.out).join(contig_map),
        binFiles(['join_flags.py'])
    )

    STAGE_SUMMARY(
        JOIN_FLAGS.out.stats.map { _cid, f -> f }.collect().map { fs -> tuple('flags', fs) },
        binFiles(['stage_summary.py'])
    )

    emit:
    smorfs = JOIN_FLAGS.out.smorfs   // tuple(chunk_id, smorfs.tsv)
}

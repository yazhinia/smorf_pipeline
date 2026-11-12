include { BAKTA; SMORFINDER; GMSC_MAPPER; GMSC_EXTRACT; CUSTOM_CALLER } from '../../modules/local/callers'
include { STAGE_SUMMARY } from '../../modules/local/report'
include { binFiles } from '../../modules/local/utils'

def chunkFastas(genomeDir) {
    return genomeDir.listFiles().findAll { f -> f.name.endsWith('.fna') }.sort { f -> f.name }
}

// Split one chunk's genomes into call batches. Batches never cross a chunk, and groupKey
// carries the batch count so each chunk is regrouped (and merged) as soon as its own
// batches finish. Small runs use one genome per task.
def batchChunk(String cid, genomeDir, int batchSize) {
    def fastas = chunkFastas(genomeDir)
    def size = fastas.size() <= params.small_run ? 1 : batchSize
    def batches = fastas.collate(size)
    def key = groupKey(cid, batches.size())
    return (0..<batches.size()).collect { i -> tuple(key, i, batches[i]) }
}

// The same for every extra caller at once: one group per chunk over all their batches, so
// the chunk is merged only when every extra caller has finished it.
def extraBatches(String cid, genomeDir, List callers) {
    def fastas = chunkFastas(genomeDir)
    def jobs = []
    callers.each { c ->
        def size = fastas.size() <= params.small_run ? 1 : ((c.batch ?: 20) as int)
        fastas.collate(size).eachWithIndex { b, i -> jobs << [c, i, b] }
    }
    def key = groupKey(cid, jobs.size())
    return jobs.collect { c, i, b -> tuple(c, key, i, b) }
}

workflow CALL_ORFS {
    take:
    genomes        // tuple(chunk_id, genomes_dir)
    chunk_fasta    // tuple(chunk_id, chunk.fna, contig_map)

    main:
    BAKTA(genomes.flatMap { cid, gdir -> batchChunk(cid, gdir, params.bakta_batch) })
    SMORFINDER(genomes.flatMap { cid, gdir -> batchChunk(cid, gdir, params.smorf_batch) })
    nonempty = chunk_fasta.filter { _cid, fna, _cmap -> fna.size() > 0 }
    GMSC_MAPPER(nonempty.map { cid, fna, _cmap -> tuple(cid, fna) })
    GMSC_EXTRACT(GMSC_MAPPER.out.join(nonempty.map { cid, _fna, cmap -> tuple(cid, cmap) }),
                 binFiles(['gmsc_extract.py']))
    CUSTOM_CALLER(genomes.flatMap { cid, gdir -> extraBatches(cid, gdir, params.extra_callers ?: []) },
                  binFiles(['run_caller.sh', 'custom_gff.py', 'pam_modules.sh']))

    bakta = BAKTA.out.groupTuple().map { key, dirs -> tuple(key.getGroupTarget(), dirs) }
    smorf = SMORFINDER.out.groupTuple().map { key, dirs -> tuple(key.getGroupTarget(), dirs) }
    // No extra callers: an empty list per chunk, so the merge's join still fires.
    extra = params.extra_callers
        ? CUSTOM_CALLER.out.groupTuple().map { key, dirs -> tuple(key.getGroupTarget(), dirs.sort { d -> d.name }) }
        : genomes.map { cid, _gdir -> tuple(cid, []) }

    all_dirs = bakta.map { _cid, dirs -> dirs }
        .mix(smorf.map { _cid, dirs -> dirs }, GMSC_EXTRACT.out.map { _cid, dir -> [dir] },
             extra.map { _cid, dirs -> dirs })
        .collect()
    STAGE_SUMMARY(all_dirs.map { ds -> tuple('call', ds) }, binFiles(['stage_summary.py']))

    emit:
    bakta = bakta               // tuple(chunk_id, [bakta dirs])
    smorf = smorf               // tuple(chunk_id, [smorfinder dirs])
    gmsc  = GMSC_EXTRACT.out    // tuple(chunk_id, gmsc dir)
    extra = extra               // tuple(chunk_id, [extra caller dirs]), [] without extra_callers
}

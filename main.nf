#!/usr/bin/env nextflow
/*
 * smorf_pipeline: genome FASTA -> candidate microproteins (small ORFs <= 100 aa),
 * using the microprotein catalog method.
 *
 *   prepare  -> plan genomes, normalize contig ids, chunk
 *   call     -> Bakta, smORFinder, GMSC-Mapper (high-quality hits), + extra_callers
 *   merge    -> consensus calls per genome, translate, keep <= max_aa
 *   combine  -> + orfipy novel intergenic ORFs, sequence hygiene
 *   flags    -> AntiFam and Popcorn flags, final tables
 *
 * --stop_after <stage> pauses after that stage; relaunch with -resume to carry on.
 */

include { PREPARE }         from './subworkflows/local/prepare'
include { CALL_ORFS }       from './subworkflows/local/call_orfs'
include { MERGE }           from './subworkflows/local/merge'
include { COMBINE }         from './subworkflows/local/combine'
include { FLAGS }           from './subworkflows/local/flags'
include { COLLECT_RESULTS } from './modules/local/report'
include { binFiles }        from './modules/local/utils'

def stages() {
    return ['prepare', 'call', 'merge', 'combine', 'flags']
}

def checkParams() {
    if (!params.input) {
        error "Missing --input: a FASTA file, a directory of FASTAs, or a samplesheet (.tsv/.csv)"
    }
    if (!(params.stop_after in stages())) {
        error "--stop_after must be one of ${stages().join(', ')}; got '${params.stop_after}'"
    }
    if (!(params.split_mode in ['file', 'record', 'regex'])) {
        error "--split_mode must be file, record or regex; got '${params.split_mode}'"
    }
    if (params.split_mode == 'regex' && !params.split_regex) {
        error "--split_mode regex needs --split_regex with one capture group"
    }
    checkExtraCallers(params.extra_callers)
    def unset = ['antifam_hmm', 'popcorn_dir'].findAll { k -> !params[k] }
    if (unset && params.stop_after == 'flags') {
        error "Not set: ${unset.join(', ')}. Copy conf/site.config.example to conf/site.config and fill in the site paths."
    }
    // PaM tools are Singularity wrappers that only see /lustre and /data/pam.
    def wd = workflow.workDir.exists() ? workflow.workDir.toRealPath().toString() : workflow.workDir.toString()
    if (workflow.profile.tokenize(',').contains('sanger') && !(wd.startsWith('/lustre/') || wd.startsWith('/data/pam/'))) {
        error "Work dir ${wd} is not on /lustre or /data/pam; the PaM tool containers cannot read it. Use -w /lustre/..."
    }
}

// Same rules as bin/list_tools.py, which explains them to the user before launch.
def checkExtraCallers(callers) {
    if (!callers) {
        return
    }
    if (!(callers instanceof List) || !callers.every { c -> c instanceof Map }) {
        error "extra_callers must be a list of {name, env, command} entries"
    }
    def reserved = ['bakta', 'smorfinder', 'gmscmapper', 'orfipy', 'pyrodigal', 'prodigal_v2.6.3']
    def names = callers.collect { c -> (c.name ?: '') as String }
    callers.each { c ->
        def name = (c.name ?: '') as String
        if (!(name ==~ /[a-z][a-z0-9_]{0,30}/) || name in reserved || names.count(name) > 1) {
            error "extra_callers: bad or repeated name '${name}' (lowercase letters, digits, _; not ${reserved.join(', ')})"
        }
        if (!(((c.env ?: 'none') as String) ==~ /module:[\w.\/,+:=-]+|(conda|sif):\/[\w.\/,+:=-]+|none/)) {
            error "extra_callers '${c.name}': env must be module:<name>, conda:/<path>, sif:/<image> or none"
        }
        if (!(c.command instanceof String) || !c.command.contains('{fasta}') || !c.command.contains('{gff}')) {
            error "extra_callers '${c.name}': command must use {fasta} and {gff}"
        }
    }
}

workflow {
    checkParams()
    def stop = stages().indexOf(params.stop_after)
    def input_path = file(params.input, checkIfExists: true).toAbsolutePath().toString()

    log.info "smorf_pipeline ${workflow.manifest.version}: input=${input_path} stop_after=${params.stop_after} outdir=${params.outdir}"

    PREPARE(channel.value(input_path))

    if (stop >= 1) {
        CALL_ORFS(PREPARE.out.genomes, PREPARE.out.chunk_fasta)

        if (stop >= 2) {
            MERGE(PREPARE.out.genomes, PREPARE.out.stats,
                  CALL_ORFS.out.bakta, CALL_ORFS.out.smorf, CALL_ORFS.out.gmsc,
                  CALL_ORFS.out.extra)

            if (stop >= 3) {
                COMBINE(PREPARE.out.genomes, MERGE.out.merged, MERGE.out.table)

                if (stop >= 4) {
                    FLAGS(PREPARE.out.genomes, COMBINE.out.smorfs, PREPARE.out.contig_map)

                    final_inputs = FLAGS.out.smorfs.map { _cid, f -> f }
                        .mix(
                            PREPARE.out.stats.map { _cid, f -> f },
                            PREPARE.out.contig_map.map { _cid, f -> f },
                            MERGE.out.table.map { _cid, _t, s -> s },
                            MERGE.out.info.map { _cid, skipped, _c -> skipped },
                            COMBINE.out.info.flatMap { _cid, dropped, cmap, _s -> [dropped, cmap] }
                        )
                        .collect()
                    COLLECT_RESULTS(final_inputs, binFiles(['collect_results.py']))
                }
            }
        }
    }
}

// Stage "call": the three ORF callers the catalog merges, plus any extra_callers.
// A genome that fails is written to status.tsv and the task still succeeds, so one bad
// genome never stops a batch (same policy as the production *_worker.sh scripts).

process BAKTA {
    tag "${chunk_key}:b${batch_idx}"
    label 'bakta'

    input:
    tuple val(chunk_key), val(batch_idx), path(fastas, stageAs: 'in/*')

    output:
    tuple val(chunk_key), path("${chunk_key}.bakta_b${batch_idx}", type: 'dir')

    script:
    // Command from bakta_worker.sh:114-124; ##FASTA strip from :147. Genomes that fail or
    // time out are retried once at the end with the longer timeout (schedule_bakta.sh:53).
    def out = "${chunk_key}.bakta_b${batch_idx}"
    """
    source ${projectDir}/bin/pam_modules.sh ${params.bakta_module}
    export OMP_NUM_THREADS=${task.cpus}
    mkdir -p ${out}
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${out}/status.tsv

    run_bakta() {   # <fasta> <timeout>; 0 on success
        local fa=\$1 limit=\$2 gid rc=0
        gid=\$(basename "\$fa" .fna)
        rm -rf "work_\$gid"; mkdir -p "work_\$gid"
        timeout "\$limit" bakta \\
            --db ${params.bakta_db} \\
            --keep-contig-headers \\
            --locus-tag "\$gid" \\
            --skip-crispr --skip-plot --skip-ori --skip-gap \\
            --threads ${task.cpus} \\
            --tmp-dir "work_\$gid" \\
            --output "work_\$gid" \\
            --prefix "\$gid" \\
            --force \\
            "\$fa" > "work_\$gid/bakta.stdout" 2> "work_\$gid/bakta.stderr" || rc=\$?
        if [ "\$rc" -eq 0 ] && [ -s "work_\$gid/\$gid.gff3" ]; then
            awk '/^##FASTA/{exit} {print}' "work_\$gid/\$gid.gff3" > "${out}/\$gid.gff3"
            n=\$(grep -cv '^#' "${out}/\$gid.gff3" || true)
            printf '%s\\tbakta\\tok\\t%s\\t%s\\n' "\$gid" "\$3" "\$n" >> ${out}/status.tsv
            rm -rf "work_\$gid"
            return 0
        fi
        [ "\$rc" -eq 0 ] && rc=no_gff3
        tail -5 "work_\$gid/bakta.stderr" >&2 || true
        echo "\$rc" > "work_\$gid/rc"
        return 1
    }

    failed=()
    for fa in in/*.fna; do
        run_bakta "\$fa" ${params.bakta_timeout} "" || failed+=("\$fa")
    done
    for fa in "\${failed[@]+"\${failed[@]}"}"; do
        gid=\$(basename "\$fa" .fna)
        echo "retrying \$gid with ${params.bakta_retry_timeout}" >&2
        if ! run_bakta "\$fa" ${params.bakta_retry_timeout} retry_ok; then
            rc=\$(cat "work_\$gid/rc" 2>/dev/null || echo unknown)
            reason="bakta_rc_\$rc"
            [ "\$rc" = "124" ] && reason=timeout
            printf '%s\\tbakta\\tfailed\\t%s\\t0\\n' "\$gid" "\$reason" >> ${out}/status.tsv
            rm -rf "work_\$gid"
        fi
    done
    """

    stub:
    def out = "${chunk_key}.bakta_b${batch_idx}"
    """
    mkdir -p ${out}
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${out}/status.tsv
    for fa in in/*.fna; do gid=\$(basename "\$fa" .fna); touch ${out}/\$gid.gff3
        printf '%s\\tbakta\\tok\\t\\t0\\n' "\$gid" >> ${out}/status.tsv; done
    """
}

process SMORFINDER {
    tag "${chunk_key}:b${batch_idx}"
    label 'smorfinder'

    input:
    tuple val(chunk_key), val(batch_idx), path(fastas, stageAs: 'in/*')

    output:
    tuple val(chunk_key), path("${chunk_key}.smorf_b${batch_idx}", type: 'dir')

    script:
    // smorfinder_worker.sh:77-146: 'smorf single', failure classes, and the empty-GFF
    // rescue for genomes with no protein <= 51 aa (smORFinder crashes instead of writing one).
    def out = "${chunk_key}.smorf_b${batch_idx}"
    """
    source ${projectDir}/bin/pam_modules.sh ${params.smorfinder_module}
    mkdir -p ${out}
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${out}/status.tsv

    for fa in in/*.fna; do
        gid=\$(basename "\$fa" .fna)
        wd="\$gid"                      # smORFinder names its outputs after the output dir
        rm -rf "\$wd"; mkdir -p "\$wd"
        log="\$gid.smorf.log"
        rc=0
        smorf single "\$fa" -o "\$wd" > "\$log" 2>&1 || rc=\$?
        if [ "\$rc" -eq 0 ] && [ -f "\$wd/\$gid.gff" ]; then
            mv "\$wd/\$gid.gff" "${out}/\$gid.gff"
            n=\$(grep -cv '^#' "${out}/\$gid.gff" || true)
            printf '%s\\tsmorfinder\\tok\\t\\t%s\\n' "\$gid" "\$n" >> ${out}/status.tsv
            rm -rf "\$wd" "\$log"
            continue
        fi
        reason=unknown
        if grep -q 'free(): invalid pointer' "\$log"; then reason=prodigal_abort
        elif grep -q 'cannot convert float infinity to integer' "\$log"; then reason=no_small_orfs_keras_progbar
        elif grep -q 'prodigal.small.faa is empty or misformatted' "\$log"; then reason=no_small_orfs
        elif grep -q 'IndexError: list index out of range' "\$log"; then reason=rename_proteins_indexerror
        elif [ "\$rc" -eq 0 ]; then reason=no_gff_produced
        fi
        n_small=NA
        [ -f "\$wd/tmp/prodigal.small.faa" ] && n_small=\$(grep -c '^>' "\$wd/tmp/prodigal.small.faa" || true)
        if [ "\$n_small" = "0" ]; then
            n_faa=\$(grep -c '^>' "\$wd/tmp/prodigal.faa" || true)
            n_ffn=\$(grep -c '^>' "\$wd/tmp/prodigal.ffn" || true)
            n_gff=\$(awk '!/^#/ && length(\$0) > 10' "\$wd/tmp/prodigal.gff" | wc -l)
            n_bare=\$(grep -c '^>\$' "\$wd/tmp/prodigal.faa" || true)
            if [ "\$n_faa" -gt 0 ] && [ "\$n_faa" = "\$n_ffn" ] && [ "\$n_faa" = "\$n_gff" ] && [ "\$n_bare" = "0" ]; then
                : > "${out}/\$gid.gff"
                printf '%s\\tsmorfinder\\tempty_ok\\tno_small_orfs\\t0\\n' "\$gid" >> ${out}/status.tsv
                rm -rf "\$wd" "\$log"
                continue
            fi
        fi
        echo "\$gid FAILED reason=\$reason rc=\$rc small_faa=\$n_small" >&2
        tail -4 "\$log" >&2 || true
        printf '%s\\tsmorfinder\\tfailed\\t%s\\t0\\n' "\$gid" "\$reason" >> ${out}/status.tsv
        rm -rf "\$wd"
    done
    """

    stub:
    def out = "${chunk_key}.smorf_b${batch_idx}"
    """
    mkdir -p ${out}
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${out}/status.tsv
    for fa in in/*.fna; do gid=\$(basename "\$fa" .fna); touch ${out}/\$gid.gff
        printf '%s\\tsmorfinder\\tok\\t\\t0\\n' "\$gid" >> ${out}/status.tsv; done
    """
}

process GMSC_MAPPER {
    tag "${chunk_id}"
    label 'gmsc'

    input:
    tuple val(chunk_id), path(chunk_fna)

    output:
    tuple val(chunk_id), path("${chunk_id}.gmsc_raw", type: 'dir')

    script:
    // 04_run_gmsc_mapper.sh:47-52. Only the three small files the high-quality filter reads
    // are kept; the ~5 GB alignment table and the rest are deleted in the task.
    def raw = "${chunk_id}.gmsc_raw"
    """
    source ${projectDir}/bin/pam_modules.sh ${params.gmsc_module}
    mkdir -p gmsc_out ${raw}
    gmsc-mapper \\
        -i ${chunk_fna} \\
        --db ${params.gmsc_db} \\
        -o gmsc_out \\
        --tool mmseqs \\
        -t ${task.cpus}
    mv gmsc_out/quality.out.smorfs.tsv gmsc_out/mapped.smorfs.faa gmsc_out/predicted.filterd.smorf.faa ${raw}/
    rm -rf gmsc_out
    """

    stub:
    """
    mkdir -p ${chunk_id}.gmsc_raw
    """
}

process CUSTOM_CALLER {
    // One extra caller from params.extra_callers. `tool` is that entry, so editing it re-runs
    // only this tool's tasks; resources come from the entry, not from base.config.
    tag "${tool.name}:${chunk_key}:b${batch_idx}"
    cpus   { (tool.cpus ?: 1) as int }
    memory { new nextflow.util.MemoryUnit((tool.memory ?: '4 GB') as String) * task.attempt }
    time   { (tool.time ?: '4h') as String }

    input:
    tuple val(tool), val(chunk_key), val(batch_idx), path(fastas, stageAs: 'in/*')
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    tuple val(chunk_key), path("${chunk_key}.x_${tool.name}_b${batch_idx}", type: 'dir')

    script:
    // The user's command goes through base64, so quotes, pipes and newlines arrive unchanged.
    def out = "${chunk_key}.x_${tool.name}_b${batch_idx}"
    def b64 = { String s -> s.bytes.encodeBase64().toString() }
    def version = tool.version_cmd ? "--version-cmd \"\$(echo ${b64(tool.version_cmd as String)} | base64 -d)\"" : ''
    """
    echo ${b64(tool.command as String)} | base64 -d > command.txt
    export SMORF_SINGULARITY_MODULE=${params.singularity_module}
    ${projectDir}/bin/run_caller.sh --name ${tool.name} --env '${tool.env ?: 'none'}' \\
        --command-file command.txt --out-dir ${out} --cpus ${task.cpus} \\
        --timeout ${tool.timeout ?: '1h'} --python ${params.python_bin} ${version} in/*.fna
    """

    stub:
    def out = "${chunk_key}.x_${tool.name}_b${batch_idx}"
    """
    mkdir -p ${out}
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${out}/status.tsv
    for fa in in/*.fna; do gid=\$(basename "\$fa" .fna); touch ${out}/\$gid.gff
        printf '%s\\t${tool.name}\\tok\\t\\t0\\n' "\$gid" >> ${out}/status.tsv; done
    """
}

process GMSC_EXTRACT {
    tag "${chunk_id}"
    label 'py_small'

    input:
    tuple val(chunk_id), path(raw), path(contig_map)
    path code, stageAs: 'code/*'     // the bin/ scripts used, so code edits re-run this step

    output:
    tuple val(chunk_id), path("${chunk_id}.gmsc", type: 'dir')

    script:
    // High-quality filter of process_gmscmapper_chunks.sh:12-17, one GFF per genome.
    """
    ${params.python_bin} ${projectDir}/bin/gmsc_extract.py \\
        --gmsc-out ${raw} \\
        --contig-map ${contig_map} \\
        --out-dir ${chunk_id}.gmsc
    """

    stub:
    """
    mkdir -p ${chunk_id}.gmsc
    printf 'genome_id\\ttool\\tstatus\\treason\\tn_records\\n' > ${chunk_id}.gmsc/status.tsv
    """
}

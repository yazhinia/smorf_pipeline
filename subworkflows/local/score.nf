// Placeholder for a later scoring stage (e.g. a translation model). Not wired in v1.
//
// Contract: take tuple(chunk_id, smorfs.tsv) from FLAGS, add score columns keyed on orf_id,
// emit tuple(chunk_id, smorfs.tsv) with the same columns plus the new ones, and write
// qc/score_summary.* through STAGE_SUMMARY like every other stage.

workflow SCORE {
    take:
    smorfs

    main:
    scored = smorfs

    emit:
    smorfs = scored
}

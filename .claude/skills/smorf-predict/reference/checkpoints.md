# Checkpoints: what to show, what is normal, what to ask

Reference genome for "normal": catalog genome gI053058 (*E. coli* isolate, 2 contigs, 5.2 Mb).
Catalog-wide figures are over 155k human-gut genomes (HRGM2).

## prepare (`qc/prepare_summary.*`)
| metric | meaning | normal / action |
|---|---|---|
| genomes_ok / genomes_failed | genomes read correctly / empty or unreadable files | failed > 0: name them; ask to fix the input or go on without them |
| genome_size_bp_min_median_max | assembly size in base pairs | bacteria 0.5-10 Mb; < 100 kb is suspicious |
| non_acgt_pct | share of bases that are N or IUPAC codes | > 1% WARN: smORFs with X are dropped later |
| genomes_over_999_contigs | very fragmented assemblies | WARN only: ORF numbering order is odd, results are fine |
| est_bakta_core_hours, gmsc_200gb_jobs | cost of the next stage | always say these before launching `call` |

The tools question (`list_tools.py` table, every tool by stage) is asked while `prepare` runs;
gene callers can be added, a newer tool for other stages is logged as a request. Ask: continue
to `call` (state the cost; the default), run the rest without stopping (state the cost of all
stages), adjust input/split, or stop.

## call (`qc/call_summary.*`)
| metric | meaning | normal / action |
|---|---|---|
| tools | the callers that ran, built-in first | extra callers listed after gmscmapper |
| bakta_ok | genomes Bakta annotated | must be all; failed genomes are skipped by the merge (reason column: timeout, bakta_rc_N) |
| smorfinder_ok | includes `empty_ok` = genome has no protein <= 51 aa | failures (prodigal_abort, ...) just mean no smORFinder support for that genome |
| gmscmapper_ok | GMSC-Mapper high-quality hits per genome | 0 hits is possible for novel taxa |
| *_records_per_genome_median | features per genome | gI053058: Bakta 5,453 lines incl. genes/RNAs; smORFinder 16; GMSC 336 |
| <extra>_ok | genomes the extra caller ran on (`ok`, or `empty_ok` = no CDS) | all failed is FAIL: wrong command or env; reasons `rc_N`, `timeout`, `no_gff`, `unknown_contigs:<id>` (GFF not on the genome's contig names) |
| <extra>_version | first line of its `version_cmd` | |

Tools table: `list_tools.py <run>` (SKILL.md 3b). Example extra caller on gI053058:
`pyrodigal_sorf` (Pyrodigal 3.7.1, `--min-gene 45`) 5,106 CDS, 703 of them <= 100 aa.

Ask: continue to `merge`, retry failed genomes (re-launch `call`; Bakta already retries once with 3 h), or stop.

## merge (`qc/merge_summary.*`)
| metric | meaning | normal |
|---|---|---|
| cds_total | all merged gene calls, any length | gI053058: 5,103 |
| microproteins_total / per_genome_median | merged calls <= 100 aa | gI053058: 706; catalog mean ~240 per genome (MAGs are smaller) |
| internal_stop | translated with a stop inside | 0 expected; > 0 is removed at `combine` |
| support:<tools> | which callers found each microprotein, % of total (extra callers appear by name) | catalog: bakta+gmscmapper 46.8%, bakta 44.5%, bakta+smorfinder 3.2%, gmscmapper only 2.8%, all three 2.1%, smorfinder only 0.5% |

Explain "support": a smORF found by more than one caller is more trustworthy.
Ask: continue to `combine`, or stop with the merged table (`chunks/*/*.merged_orf_table.tsv` when stopped here).

## combine (`qc/combine_summary.*`, `qc/combine_aa_length.png`)
| metric | meaning | normal |
|---|---|---|
| smorfs_merged_catalog | merged smORFs that passed hygiene | gI053058: 706 |
| smorfs_orfipy_novel | every 15-300 nt ORF in unannotated DNA | gI053058: ~16,600; this is a discovery pool, mostly not real proteins |
| dropped_hygiene | removed for X, internal stop or odd letters | ~0.1% |
| partial:<source> | smORFs with no start or stop codon, or starting at a contig end; the real protein may be longer (a fragment of a big protein can pass as ≤100 aa) | 0 for a complete genome (gI053058: 0); more in fragmented drafts and MAGs; > 1% WARN |
| aa_length_median:<source> | median protein length | |

Figure: x = protein length (aa); y = % of that source's smORFs at each length; one line per source.
Ask: continue to `flags`, or stop.

## flags (`qc/flags_summary.*`, `qc/flags_popcorn.png`)
| metric | meaning | normal |
|---|---|---|
| antifam_hit_pct:<source> | % matching a known spurious-ORF family (AntiFam) | low single digits; > 25% WARN |
| popcorn_coding_pct:<source> | % Popcorn calls coding | higher for merged_catalog than orfipy_novel |
| popcorn_missing_pct:<source> | Popcorn skipped the sequence | left empty, not "non-coding" |

Figure: x = Popcorn coding probability (0.05 bins); y = % of that source's smORFs per bin.

## Final outputs (`results/`)
| file | content |
|---|---|
| smorfs.tsv | one row per smORF: ids, position, length, `partial` (`5prime`/`3prime`/`both`/`none`), callers (`tools`, `support_n`), `source`, AntiFam and Popcorn flags, protein sequence, sample name, original contig name |
| smorfs.faa / smorfs.gff3 | same smORFs as FASTA, and as GFF3 on the user's own contig names |
| dropped.tsv | smORFs removed by hygiene, with the reason |
| per_genome_summary.tsv | one row per genome: input stats, merge status, smORF counts |
| genome_map.tsv / contig_map.tsv | pipeline ids back to the user's names |
| chunks/<chunk>/ | per-chunk tables (the only copy above `concat_max_genomes` genomes) |
| qc/ | checkpoint summaries and figures |
| run_log/ | Nextflow report, timeline and trace; plan.tsv (genome ids to input files) |

Outside `results/`, `<run>/run_log/` holds `decisions.tsv` (every question and answer of the run:
time, stage, question, choice, the user's comment, params change, then options offered, user,
run folder, params and code version) and the Nextflow head-job logs `head_<stage>_<time>.log`.

## What an option change re-runs
| option | first stage that re-runs |
|---|---|
| input, split_mode, split_regex, id_prefix, chunk_size | prepare (everything) |
| bakta_batch, smorf_batch | call: Bakta / smORFinder for every genome (costly; avoid) |
| extra_callers: add or edit one entry | call (that tool only), then merge onward |
| extra_callers: remove an entry | merge onward |
| bakta_*, smorfinder_module, gmsc_* | call (that tool) and everything after |
| max_aa | merge |
| orf_min, orf_max, orf_starts, orf_stops, assets/hygiene.yaml | combine |
| antifam_*, hmmer_module, popcorn_* | flags |
| max_chunks | nothing; extra chunks are added |
| publish_*, concat_max_genomes, stop_after | nothing (output only) |
| any code edit in bin/*.py | only the steps using that script and what follows them |
| edit to plan_inputs.py / normalize_chunk.py | everything (new genome files), including Bakta and GMSC-Mapper |

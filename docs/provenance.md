# Where each step comes from

Paths relative to the microprotein catalog project folder.
`S` = `orf_prediction_run/submission_scripts`, `K` = `.codex/skills/hrgm2-orf-workflow/scripts`.

| Pipeline step | Ported from | Change |
|---|---|---|
| `bin/plan_inputs.py`, `bin/normalize_chunk.py` | `K/00_prepare_inputs.py` | no HRGM2 metadata; `gU` ids; record/regex split modes |
| `BAKTA` | `S/bakta_scheduler/bakta_worker.sh:114-147`, `schedule_bakta.sh` (4 threads, 3 h retry) | in-task retry; module 1.11.0 (1.11.4 retired) |
| `SMORFINDER` | `S/smorfinder_scheduler/smorfinder_worker.sh:77-146` | none |
| `GMSC_MAPPER`, `bin/gmsc_extract.py` | `K/04_run_gmsc_mapper.sh:47-52`, `K/process_gmscmapper_chunks.sh:12-17` | exact id match instead of `grep -F`; genome from contig map |
| `bin/merge_gffs.py` | `K/merge_gffs.py` | none (byte-identical) |
| `bin/merge_chunk.py` | `S/merge_tools.py:129-246` | per-genome inputs instead of batch files |
| `ORFIPY_NOVEL`, `bin/filter_novel_orfs.py` | `S/orfipy_scheduler/orfipy_worker.sh:30-34,123-129`, `filter_novel_orfs.py` | import path only |
| `bin/combine_chunk.py`, `bin/hygiene.py`, `assets/hygiene.yaml` | `S/cluster_scheduler/combine_batch_orfs.py`, `S/classifier_dataset/rmambiguous_seq.py:32-51`, `lf_config.yaml` hygiene block | genomes from merge summary; `partial` flag is new |
| `bin/extract_seqs.py` | `S/cluster_scheduler/extract_rep_nucleotide.py:93-98` | all distinct sequences, not reps |
| `ANTIFAM` | `S/cluster_scheduler/antifam_worker.sh:98-114` | `-Z` fixed at catalog size |
| `POPCORN` | `S/cluster_scheduler/run_popcorn.sh:56-59` | env python called directly |
| `bin/join_flags.py` | `S/cluster_scheduler/collect_antifam.py:93-121`, `join_popcorn_labels.py` | per smORF |

Verified: `tests/compare_to_catalog.py fast-merge --genome gI053058` reproduces the catalog's
706 merged microprotein rows exactly (all 11 columns).

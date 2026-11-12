---
name: smorf-predict
description: Predict microproteins (small ORFs, <=100 aa) from bacterial genome FASTA files with the smorf_pipeline (catalog method - Bakta, smORFinder, GMSC-Mapper merge, orfipy novel ORFs, AntiFam and Popcorn flags) on the Sanger farm, step by step with the user - inspect the input, agree options, run one stage at a time, review each checkpoint summary, and log every decision for the developers. Shows every tool by stage (gene callers, intergenic smORFs, flags) with versions, and can add the user's own gene caller (install, test on one genome, merge). Use when someone wants smORFs/microproteins called from one genome, a multi-FASTA, or a folder of assemblies, wants to add or swap a gene caller, or wants to continue, adjust or check such a run.
---

# smORF prediction, one checkpoint at a time

You run the pipeline **with** the user. Each stage ends at a checkpoint: show the summary in
plain words, ask what to do next, log the answer, then go on. Never run two stages without
the user's say-so (choosing "Run the rest without stopping", step 5c, is that say-so), and
never change an option silently. While a stage runs, use the wait (step 5a).

The user is usually in their own working folder, not in the pipeline. Before anything else,
find `PIPE`, the pipeline copy (the folder holding `main.nf`); scripts are in `$PIPE/bin`:
1. `realpath "${CLAUDE_SKILL_DIR}/../../.."`: this skill sits at `<PIPE>/.claude/skills/smorf-predict`,
   and `realpath` also follows the `~/.claude/skills/smorf-predict` link the README sets up.
2. If that folder has no `main.nf`, use `$SMORF_PIPELINE` if it is set.
3. Otherwise ask the user where the pipeline copy is.
Check `$PIPE/main.nf` exists and use full `$PIPE/...` paths from then on.

Read `reference/checkpoints.md` (next to this file) before the first checkpoint; it says what
each summary metric means, what is normal, and what to ask.

## Hard rules
- **Nothing runs on the head node.** Inspection goes through `bsub -I -q normal -n 1 -M 2000 -R 'select[mem>2000] rusage[mem=2000]' ...`;
  the pipeline itself goes through `$PIPE/bin/smorf_launch.sh`, which submits the Nextflow head job.
- The run folder must be on `/lustre` or `/data/pam` (PaM tool containers cannot see `/tmp`, `/nfs` home or other paths).
- Never write into `$PIPE`; everything a run makes goes in its run folder, and tools the user
  adds go in a shared tools folder (default `<working folder>/smorf_tools`).
- Every user decision is logged with `$PIPE/bin/log_decision.py` (step 6). This log is how the
  pipeline gets improved, so include the user's own words in `--comment` when they give a reason.
- Explain results in simple words; for any table or figure, say what the rows, columns or axes are.

## 1. Input and run folder
Ask for (unless already given): the input path, and a run folder name. Default run folder:
`<current working folder>/smorf_runs/<name>` when the working folder is on `/lustre` or
`/data/pam`; otherwise ask for one there. One run folder per genome set; a new set gets a
new folder.

Several genome sets can run in parallel: give each its own run folder and launch them one
after another; they share nothing but `$PIPE`. Keep a short list of the user's runs (folder,
input, current stage) and check each with `smorf_status.sh`. Never launch two runs in the same
run folder at once.

Input can be one FASTA (`.fa/.fna/.fasta`, optionally `.gz`), a directory of FASTAs (one
genome per file), or a samplesheet (`.tsv`/`.csv` with columns `sample`, `fasta`, optional `genome_id`).

## 2. Inspect
```bash
bsub -I -q normal -n 1 -M 2000 -R 'select[mem>2000] rusage[mem=2000]' python3 $PIPE/bin/inspect_inputs.py <input>
```
Tell the user what it found (files, records, sizes) and its suggested `split_mode`; keep the
example genome and contig names for the id question in step 3.
- One record, or a folder/samplesheet → `split_mode: file`; no question needed.
- A multi-record FASTA → **ask**: "Is this one genome in several contigs, one genome per record,
  or several genomes whose record names share a prefix?" → `file`, `record` or `regex`
  (with the suggested `split_regex`, shown to the user).
- **Ask** whether the genomes are bacterial or archaeal. The method uses genetic code 11 and
  bacterial gene callers; for anything else, say it is not supported and stop.
- Total under 100 kb, or many non-ACGT letters: point it out before going on.

## 3. Options → params.yaml
Write `<run>/params.yaml` from `reference/params_template.yaml`. Then explain the options that
matter, each in two or three plain sentences with an example from the user's own input, and ask
them together in one AskUserQuestion (recommended option first). Show the final `params.yaml`.

**Genome ids (`id_prefix`).** Explain: the pipeline gives every genome a short id, its contigs
and smORFs follow from it, and the original names are never lost. Say why it helps: genome and
contig names are often long or contain dots, spaces or `|` (e.g. `GCF_000005845.2_ASM584v2_genomic`,
`NODE_12_length_52341_cov_12.5`), which some tools cut short or reject, and short ids keep every
tool's output in step. If the inspect output flagged long or odd names, say so with one of them.
Show a small table from the inspect output's names:

| your name | pipeline id |
|---|---|
| `GENOME148803` (genome) | `demo000001` |
| first contig, second contig | `demo000001_c001`, `demo000001_c002` |
| a smORF on it | `demo000001_ORF00042` |

Then: "`genome_map.tsv` and `contig_map.tsv` link ids to your names, and `smorfs.tsv` also has
your `sample_name` and `original_contig` on every row."
Options: (1) a short prefix for this genome set, suggested from the folder or project name (1-10
letters/digits, starting with a letter) (recommended); (2) the default `gU`; (3) their own id
per genome, via a samplesheet with a `genome_id` column (3-16 letters/digits, starting with a
letter, unique).

**Pilot run (`max_chunks`), only over ~200 genomes.** "A pilot runs the first 500 genomes all
the way through, so you can check results and cost before committing to the whole set; the
full run then reuses everything the pilot finished." Example from their count: 12,000 genomes →
pilot of 500, then the remaining 11,500. Options: pilot first (`max_chunks: 1`) (recommended),
full run.

**Per-tool files (`publish_intermediates`).** "Keep each caller's own GFF (Bakta, smORFinder,
GMSC-Mapper, orfipy, extra callers) in `results/intermediates/`: useful to compare callers or
check one call, but about 4 MB per genome." Example from their count: 10,000 genomes → ~40 GB.
Options: no (recommended above ~100 genomes), yes.

**Share decisions (`feedback_dir`).** "Your answers to these questions (not your genomes or
results) can also go to a shared log the developers read to improve the pipeline."
Options: no (default), yes (then ask for or confirm the shared folder path).

## 3b. Tools
While `prepare` runs (5a), before any costly stage (and whenever the user asks which tools
run), show every prediction tool the run will use:
```bash
bsub -I -q normal -n 1 -M 2000 -R 'select[mem>2000] rusage[mem=2000]' $PY $PIPE/bin/list_tools.py <run>
```
`$PY` is `python_bin` from `$PIPE/conf/site.config`, else `$PIPE/nextflow.config` (it has PyYAML;
the system `python3` does not):
`PY=$(sed -n "/^ *python_bin *=/{s/^ *python_bin *= *'\([^']*\)'.*/\1/p;q}" $(ls $PIPE/conf/site.config 2>/dev/null) $PIPE/nextflow.config)`.
Rows are tools, grouped by `category` (and the `stage` that runs them):
- **smORF prediction** (`call`): Bakta, smORFinder, GMSC-Mapper and any extra gene caller;
  their merged calls are the `merged_catalog` smORFs.
- **intergenic smORF** (`combine`): orfipy, ORFs in DNA with no annotated gene (`orfipy_novel`).
- **coding potential** (`flags`): Popcorn, a coding probability per smORF.
- **spurious tag** (`flags`): AntiFam, known spurious-ORF families.
Coding potential and spurious tag label smORFs but remove none.

Columns: `type` (built-in or extra), `software` (version: module, conda env, image or git),
`database`, `runs` (per genome or per chunk), `resources`, `finds` (what it calls). Explain it
in one sentence per stage, then ask **while `prepare` runs** (step 5a), as one question:
"Use these tools (default), or bring in a newer or better tool?"
- **Default:** the tools shown. In the option's description, list them by category, e.g.
  "smORF prediction: Bakta, smORFinder, GMSC-Mapper · intergenic smORF: orfipy · coding
  potential: Popcorn · spurious tag: AntiFam". If the user picks the default, gives no answer, or does not
  mention tools, go on with them; never add or change a tool unasked.
- **A newer gene caller** (stage `call`): follow "Adding a gene caller" below before `call`.
- **A newer tool for `combine` or `flags`**: these cannot be swapped by the user yet. Say so,
  log the tool they named and why with `--stage feedback` (their words in `--comment`), and go
  on with the default.

## 4. Run a stage
```bash
$PIPE/bin/smorf_launch.sh <run> <stage>          # stage: prepare | call | merge | combine | flags
```
Order: `prepare` → `call` → `merge` → `combine` → `flags`. Launch the next stage only after the
user agrees at the checkpoint, or launch `flags` once they chose the no-stop run (5c).
Relaunching resumes automatically, so finished tasks are reused.

Before `call` (tools already shown during `prepare`, 3b), give the cost from `qc/prepare_summary.tsv` (`est_bakta_core_hours`, `gmsc_200gb_jobs`)
and the rough wall time: Bakta ~20 min per genome (4 at a time per job, many jobs in parallel);
GMSC-Mapper one 16-core, 200 GB job per 500-genome chunk (normal queue), even for one genome.

## 5. Wait, then show the checkpoint
### 5a. Use the wait
Right after a launch, start a background watch on `$PIPE/bin/smorf_status.sh <run>` (exit 3 =
running, 0 = done, 1 = failed): a background Bash until-loop that exits when the code is not 3;
never poll in the foreground. Tell the user roughly how long it takes (prepare: minutes; call:
~20-40 min for a few genomes, hours for thousands; merge, combine, flags: minutes to an hour),
and that they need not watch: they can keep talking, start another genome set in its own run
folder, or close Claude (the job keeps running; later "check my run <name>" picks it up).
If the head job shows `PEND`, the farm is busy: say the run is waiting for a free slot (the work
itself takes seconds to minutes for a small run). Runs of <= 50 genomes already use the
`normal` queue and run the light steps inside the head job to cut these waits.
Then use the wait for what does not need this stage's result:
- during `prepare`: show the tools table and ask the tools question (3b). A newer gene caller
  can be found or installed and its entry written now; its one-genome test waits for `prepare`.
- during `call`: say in two lines what the next checkpoints will show; answer questions.
Never ask the checkpoint question before its summary exists. If the stage ends while the user
is still answering, finish that question first, then show the checkpoint.

### 5b. Show the checkpoint
When done, read `<run>/results/qc/<stage>_summary.md` (and `_summary.tsv` for detail).
At the `call` checkpoint, start with the `tools` line: which callers ran, and the version each
extra caller reported.
Show it as a short table, explain each line in one plain sentence, and call out every WARN/FAIL
with its likely cause (see `reference/checkpoints.md`). Show the stage's PNG if there is one.
On failure: read the failed task's `.command.err` (path in the status output) and explain it.

Then ask (AskUserQuestion, at most 4 options; "Other" is always there), recommended first:
- **Continue** to the next stage
- **Run the rest without stopping** (5c; not at the last stage)
- **Retry / investigate** when something is WARN/FAIL (e.g. failed genomes), otherwise
  **Adjust** an option (then go to step 7)
- **Stop here** (results so far stay in `<run>/results`)
At `prepare`, if the tools question is still unanswered (e.g. the user came back later), ask it
first (3b), then this question.

### 5c. Run the rest without stopping
For users who do not want to wait at each checkpoint. Before launching, say plainly: every
remaining stage runs even if an early one looks wrong, so it suits small runs or data like a
run that already worked; give the remaining cost (as before `call`). A gene caller added during
the wait gets its one-genome test first (Adding a gene caller, step 4). Log the choice, then:
```bash
$PIPE/bin/smorf_launch.sh <run> flags
```
Watch it with Monitor: one event for each new `<stage>_summary.md` (newer than
`<run>/run_state.json`) and one when the run ends; re-arm when the monitor expires (call can
take hours). For each summary, post one line: stage, OK/WARN/FAIL, its key number. If a stage
is FAIL, tell the user at once and ask whether to stop the run (`bkill -J smorf_flags_<run name>`)
or let it finish. At the end, review each remaining checkpoint in order (5b, short, WARN/FAIL
first, figures shown), without asking to continue, then go to step 8.

## 6. Log the decision
```bash
python3 $PIPE/bin/log_decision.py --run-dir <run> --stage <stage> \
  --question "<what you asked>" --options "<the choices>" --choice "<their answer>" \
  [--comment "<their reason or remark, verbatim>"] [--params-diff "<key: old -> new>"]
```
Log the input/options decisions of steps 2-3 the same way with `--stage setup`.
`<run>/run_log/decisions.tsv` is the run's one log. When the user asks what was decided, show
it as a table of its first six columns (time, stage, question, choice, comment, params change).

## 7. Adjusting
Edit `params.yaml`, show the diff, and say which stages will re-run (table in
`reference/checkpoints.md`, "What an option change re-runs"). Relaunch from the earliest affected
stage's checkpoint target (the stage the user was at). Log it with `--params-diff`.

## Adding a gene caller
A new caller's CDS calls join the merge after Bakta, GMSC-Mapper and smORFinder: it adds new
smORFs and raises `support_n` where it agrees. Say plainly that results then differ from the
published catalog method. Steps, each confirmed with the user and logged (`--stage call`):
1. **Which tool**, and its command line. It must read a genome FASTA and write a GFF with `CDS`
   rows (gene finders like Prodigal/Pyrodigal do); other outputs need a conversion in the
   command. Pick a caller `name` (lowercase letters, digits, `_`; not a built-in name).
2. **Find or install it.** A PaM module (`module avail <tool>` in a `bsub -I` job) needs no
   install: `env: module:<name/version>`. Otherwise install into the shared tools folder:
   ```bash
   bsub -I -q normal -n 2 -M 8000 -R 'select[mem>8000] rusage[mem=8000]' \
     $PIPE/bin/install_tool.sh <tools_dir> <folder> conda <bioconda package[=version]>
   #   ... pip <pypi package>     or     ... sif docker://<image>
   ```
   It prints the `env:` value and records versions in `<tools_dir>/<folder>/packages.txt`.
   Check the tool's options with `--help` (bsub -I) before writing the command.
3. **Write the entry** in `params.yaml` (`extra_callers:`, see the params template), show it, and
   check it with `list_tools.py` (3b), which explains any problem. Placeholders: `{fasta}`,
   `{gff}` (required), `{cpus}`, `{gid}`, `{tmp}`; never quote them.
4. **Test on one genome** (needs `prepare` done), the same script the pipeline runs:
   ```bash
   bsub -I -q normal -n 1 -M 4000 -R 'select[mem>4000] rusage[mem=4000]' $PIPE/bin/test_caller.sh <run> <name>
   ```
   Show the status, CDS count, how many are <= 100 aa, and the first rows. A `failed` status or
   0 CDS: read the tool log lines it prints, fix the command, test again.
5. Launch (or relaunch) `call`. The built-in callers' finished tasks are reused; only the new
   tool runs, then `merge` onward re-runs. Removing a caller re-runs `merge` onward only.

## 8. Finish
After `flags`, show the final files in `<run>/results/` and what each one holds
(`reference/checkpoints.md`, "Final outputs"). Point out the main caveats:
- `orfipy_novel` smORFs are an unlabelled discovery pool, not predictions of real proteins.
- AntiFam and Popcorn are flags, not filters; there is no translation score yet.
- `partial` rows (not `none`) may be fragments of longer proteins cut by a contig end;
  suggest filtering them out for analyses that need complete proteins.

Then ask one last question: "Anything the pipeline should do differently?" and log the answer
with `--stage feedback`.

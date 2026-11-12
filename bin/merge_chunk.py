#!/usr/bin/env python3
"""Merge Bakta + GMSC-Mapper + smORFinder (+ extra callers) for every genome in one chunk.

    merge_chunk.py --genome-stats chunk.genome_stats.tsv --genomes-dir chunk.genomes \
        --bakta-dirs b0 b1 .. --smorf-dirs s0 s1 .. --gmsc-dir gmsc --prefix chunk_000001 \
        [--extra-tools name1 name2 --extra-dirs <chunk>.x_name1_b0 ..]

Port of merge_tools.py:198-246 with the batch/file-splitting plumbing removed (each tool
already gives one GFF per genome here). Per genome:
  - Bakta is required and goes first: merge_gff_files takes the genome id from the first
    file (merge_tools.py:220-227); genomes without a Bakta GFF are skipped and listed.
  - input order is Bakta, GMSC-Mapper, smORFinder, the order the catalog merged in, which
    decides which record survives an exact-duplicate tie (merge_gffs.py rule 1). Extra
    callers (params extra_callers, already cleaned by custom_gff.py) come last, in the order
    given, so without them the merge is the catalog's.
  - smORFinder lines are restricted to '<genome>_c..' feature lines, as merge_tools.py's
    split_gff_by_genome did (:83-105).
  - ORFs are translated with table 11, trailing stop removed, and kept if <= MAX_AA
    (merge_tools.py:43, :154-161).

Writes <prefix>.merged/<genome>.gff, <prefix>.merged_orf_table.tsv,
<prefix>.merged_summary.tsv, <prefix>.merge_skipped.tsv and <prefix>.merge_combos.tsv
(microproteins per tool combination, read by the stage summary).
"""
import argparse
import re
import sys
from collections import Counter
from pathlib import Path

from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_gffs import merge_gff_files  # noqa: E402

TABLE_HEADER = ["genome_id", "orf_id", "contig", "start", "end", "strand",
                "aa_length", "sequence", "tools", "support_n", "label"]


def load_contigs(fna):
    """contig id (c001) -> sequence; same keying as merge_tools.py:129-143."""
    seqs, name, parts = {}, None, []
    with open(fna) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0].rsplit("_", 1)[-1]
                parts = []
            else:
                parts.append(line.strip())
    if name:
        seqs[name] = "".join(parts)
    return seqs


def translate(contig_seq, start, end, strand):
    """Translate a 1-based inclusive CDS interval; returns protein without stop."""
    nt = Seq(contig_seq[start - 1:end])
    if strand == "-":
        nt = nt.reverse_complement()
    if len(nt) % 3:
        nt = nt[:len(nt) - len(nt) % 3]
    return str(nt.translate(table=11)).rstrip("*")


def find(dirs, name):
    for d in dirs:
        p = Path(d) / name
        if p.exists():
            return p
    return None


def genome_feature_lines(src, dst, gid):
    """Keep only '<gid>_c..' feature lines, as merge_tools.split_gff_by_genome did."""
    prefix = f"{gid}_c"
    with open(src) as fin, open(dst, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith(prefix):
                fout.write(line if line.endswith("\n") else line + "\n")


def ok_genomes(stats_path):
    with open(stats_path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["status"] == "ok":
                yield row["genome_id"]


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genome-stats", required=True)
    ap.add_argument("--genomes-dir", required=True, type=Path)
    ap.add_argument("--bakta-dirs", nargs="+", required=True)
    ap.add_argument("--smorf-dirs", nargs="*", default=[])
    ap.add_argument("--gmsc-dir")
    ap.add_argument("--extra-tools", nargs="*", default=[])
    ap.add_argument("--extra-dirs", nargs="*", default=[])
    ap.add_argument("--max-aa", type=int, default=100)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    merged_dir = Path(f"{args.prefix}.merged")
    tmp_dir = Path(f"{args.prefix}.merge_tmp")
    merged_dir.mkdir(exist_ok=True)
    tmp_dir.mkdir(exist_ok=True)
    gmsc_dirs = [args.gmsc_dir] if args.gmsc_dir else []
    # CUSTOM_CALLER names its dirs <chunk>.x_<tool>_b<batch>
    extra_dirs = {t: [d for d in args.extra_dirs
                      if re.fullmatch(rf"x_{re.escape(t)}_b\d+", Path(d).name.split(".", 1)[-1])]
                  for t in args.extra_tools}
    n_tools = 3 + len(args.extra_tools)
    support_cols = [f"n_{k}tool" for k in range(1, n_tools + 1)]

    skipped = []
    n_genomes = 0
    combos = Counter()     # tool combination of each microprotein, for the stage summary
    with open(f"{args.prefix}.merged_summary.tsv", "w") as summary, \
            open(f"{args.prefix}.merged_orf_table.tsv", "w") as table:
        summary.write("\t".join(["genome_id", "n_cds", *support_cols,
                                 "n_microprotein", "n_internal_stop"]) + "\n")
        table.write("\t".join(TABLE_HEADER) + "\n")

        for gid in ok_genomes(args.genome_stats):
            bakta = find(args.bakta_dirs, f"{gid}.gff3")
            if bakta is None:
                skipped.append((gid, "no_bakta_gff"))
                continue
            inputs = [str(bakta)]
            gmsc = find(gmsc_dirs, f"{gid}.gff")
            if gmsc is not None:
                inputs.append(str(gmsc))
            smorf = find(args.smorf_dirs, f"{gid}.gff")
            if smorf is not None:
                smorf_clean = tmp_dir / f"{gid}.gff"
                genome_feature_lines(smorf, smorf_clean, gid)
                inputs.append(str(smorf_clean))
            for tool in args.extra_tools:
                extra = find(extra_dirs[tool], f"{gid}.gff")
                if extra is not None:
                    inputs.append(str(extra))

            kept = merge_gff_files(inputs, str(merged_dir / f"{gid}.gff"))
            n_genomes += 1

            contigs = load_contigs(args.genomes_dir / f"{gid}.fna")
            support = {k: 0 for k in range(1, n_tools + 1)}
            n_micro = n_stop = 0
            for orf_id, contig, start, end, strand, tools in kept:
                support[len(tools)] = support.get(len(tools), 0) + 1
                protein = translate(contigs[contig], start, end, strand)
                if len(protein) > args.max_aa:
                    continue
                n_micro += 1
                n_stop += "*" in protein
                combos[",".join(tools)] += 1
                table.write("\t".join([
                    gid, orf_id, contig, str(start), str(end), strand,
                    str(len(protein)), protein, ",".join(tools), str(len(tools)), "",
                ]) + "\n")
            summary.write("\t".join(map(str, [gid, len(kept), *support.values(),
                                              n_micro, n_stop])) + "\n")

    with open(f"{args.prefix}.merge_skipped.tsv", "w") as fh:
        fh.write("genome_id\treason\n")
        for gid, reason in skipped:
            fh.write(f"{gid}\t{reason}\n")
    with open(f"{args.prefix}.merge_combos.tsv", "w") as fh:
        fh.write("tools\tn_microprotein\n")
        for tools, n in combos.most_common():
            fh.write(f"{tools}\t{n}\n")
    for p in tmp_dir.iterdir():
        p.unlink()
    tmp_dir.rmdir()
    print(f"{args.prefix}: merged {n_genomes}, skipped {len(skipped)}")


if __name__ == "__main__":
    main()

#!/usr/bin/env python3
"""Join AntiFam and Popcorn flags onto every smORF of one chunk.

    join_flags.py --metadata m.tsv --faa f.faa --nt-map n.tsv --antifam a.tbl \
        --popcorn p.tsv --contig-map c.tsv --prefix chunk_000001

AntiFam: best hit (highest bit score) per protein, as collect_antifam.py:93-121.
  A table without hmmsearch's '# [ok]' trailer is truncated, and the step fails.
  antifam_hit = 1 means "matches a known spurious-ORF family"; it is a flag, not a filter.
Popcorn: coding probability and call per nucleotide sequence (join_popcorn_labels.py).
  Popcorn skips some inputs by its own rules; those stay empty, never "non-coding".

Writes <prefix>.smorfs.tsv: the combine metadata columns, the flags, the protein
sequence, the sample name and the original contig name from the user's FASTA; and
<prefix>.flags_stats.tsv: counts per source and a Popcorn probability histogram.
"""
import argparse
import sys
from collections import Counter


def read_tsv(path):
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        for line in f:
            yield dict(zip(header, line.rstrip("\n").split("\t")))


def load_antifam(path):
    with open(path, "rb") as fh:
        tail = fh.read()[-4096:].decode("utf-8", "replace").rstrip("\n")
    if not tail.endswith("# [ok]"):
        sys.exit(f"{path}: no '# [ok]' trailer, hmmsearch output is incomplete")
    best, n_models = {}, Counter()
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            f = line.split(None, 6)
            if len(f) < 6:
                continue
            target, name, acc, evalue, score = f[0], f[2], f[3], f[4], f[5]
            n_models[target] += 1
            if target not in best or float(score) > best[target][0]:
                best[target] = (float(score), acc, name, score, evalue)
    return best, n_models


def load_popcorn(path):
    calls = {}
    for row in read_tsv(path):
        calls[row["Sequence ID"].split()[0]] = (row["Coding Probability"], row["Prediction"])
    return calls


def load_sequences(path):
    seqs, name = {}, None
    with open(path) as fh:
        for line in fh:
            if line.startswith(">"):
                name = line[1:].split()[0]
            else:
                seqs[name] = line.strip()
    return seqs


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    for a in ("--metadata", "--faa", "--nt-map", "--antifam", "--popcorn", "--contig-map",
              "--prefix"):
        ap.add_argument(a, required=True)
    args = ap.parse_args()

    best, n_models = load_antifam(args.antifam)
    calls = load_popcorn(args.popcorn)
    nt_of = {r["orf_id"]: r["nt_id"] for r in read_tsv(args.nt_map)}
    seqs = load_sequences(args.faa)
    contig_info = {}
    for r in read_tsv(args.contig_map):
        contig_info[r["normalized_contig_id"]] = (r["sample_name"], r["original_contig_id"])

    stats = Counter()      # (kind, source, value) -> n, read by the stage summary
    cols_out = ["genome_id", "orf_id", "source", "contig", "start", "end", "strand",
                "aa_length", "partial", "tools", "support_n", "seq_sha1",
                "antifam_hit", "antifam_n_models", "antifam_acc", "antifam_name",
                "antifam_bitscore", "antifam_evalue",
                "popcorn_probability", "popcorn_prediction",
                "sequence", "sample_name", "original_contig"]
    with open(f"{args.prefix}.smorfs.tsv", "w") as out:
        out.write("\t".join(cols_out) + "\n")
        for r in read_tsv(args.metadata):
            hit = best.get(r["seq_sha1"])
            af = (["1", str(n_models[r["seq_sha1"]]), hit[1], hit[2], hit[3], hit[4]]
                  if hit else ["0", "0", "", "", "", ""])
            pc = calls.get(nt_of[r["orf_id"]])
            sample, orig = contig_info[f"{r['genome_id']}_{r['contig']}"]
            out.write("\t".join(
                [r[c] for c in cols_out[:12]] + af + (list(pc) if pc else ["", ""])
                + [seqs[r["orf_id"]], sample, orig]) + "\n")
            src = r["source"]
            stats[("count", src, "smorfs")] += 1
            stats[("count", src, "antifam_hit")] += bool(hit)
            if pc:
                stats[("count", src, "popcorn_coding")] += pc[1].lower().startswith("cod")
                bin_lo = min(int(float(pc[0]) * 20), 19) / 20
                stats[("popcorn_bin", src, f"{bin_lo:.2f}")] += 1
            else:
                stats[("count", src, "popcorn_missing")] += 1

    with open(f"{args.prefix}.flags_stats.tsv", "w") as s:
        s.write("kind\tkey\tvalue\tn\n")
        for (kind, src, value), k in sorted(stats.items()):
            s.write(f"{kind}\t{src}\t{value}\t{k}\n")
    totals = Counter()
    for (kind, _src, value), k in stats.items():
        if kind == "count":
            totals[value] += k
    print(f"{args.prefix}: " + " ".join(f"{k}={v}" for k, v in sorted(totals.items())))


if __name__ == "__main__":
    main()

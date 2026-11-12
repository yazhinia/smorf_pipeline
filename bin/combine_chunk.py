#!/usr/bin/env python3
"""Combine merged-catalog and orfipy novel-intergenic microproteins for one chunk.

    combine_chunk.py --merged-table t.tsv --merged-summary s.tsv --orfipy-dir orfipy \
        --genomes-dir genomes --hygiene hygiene.yaml --prefix chunk_000001

Port of combine_batch_orfs.py:150-262. Per genome, in merge order:
  - merged-catalog rows are copied through (source=merged_catalog)
  - orfipy novel rows are translated from the genome FASTA (source=orfipy_novel)
  - every sequence passes the hygiene gate first (hygiene.filterseq); failures go to
    <prefix>.dropped.tsv with the reason
  - every row gets a `partial` flag (partial_ends): new, the catalog kept such ORFs unmarked

Genomes come from the merge summary, so a genome with no merged microprotein still gets
its orfipy rows (the catalog took genomes from the merged table rows instead).

Writes <prefix>.smorfs.faa, <prefix>.smorfs_metadata.tsv, <prefix>.dropped.tsv,
<prefix>.combine_map.tsv and <prefix>.combine_stats.tsv (length histogram and drop
reasons and partial counts, read by the stage summary).
"""
import argparse
import hashlib
import sys
from collections import Counter, defaultdict
from pathlib import Path

import yaml
from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from hygiene import filterseq  # noqa: E402
from merge_chunk import load_contigs  # noqa: E402

SHA1_CHARS = 20  # 80 bits -- collision-free at 4e8 distinct sequences
START_CODONS = {"ATG", "GTG", "TTG"}
STOP_CODONS = {"TAA", "TAG", "TGA"}   # table 11
_COMP = str.maketrans("ACGT", "TGCA")


def translate_nt(contig_seq, start, end, strand):
    """Same convention as merge_chunk.translate (combine_batch_orfs.py:80-93)."""
    nt = Seq(contig_seq[start - 1:end])
    if strand == "-":
        nt = nt.reverse_complement()
    if len(nt) % 3:
        nt = nt[:len(nt) - len(nt) % 3]
    return str(nt.translate(table=11)).rstrip("*")


def partial_ends(contig_seq, start, end, strand):
    """Which end of an ORF may be cut off: none, 5prime, 3prime or both.

    The 5' end is open without a start codon, or when it lies in the contig's first or last
    codon (the gene may begin off the assembly); the 3' end is open without a stop codon.
    Gene callers report genes running off contig ends this way.
    """
    first = contig_seq[start - 1:start + 2].upper()
    last = contig_seq[end - 3:end].upper()
    at_edge = start <= 3 if strand == "+" else end > len(contig_seq) - 3
    if strand == "-":
        first, last = last.translate(_COMP)[::-1], first.translate(_COMP)[::-1]
    five = at_edge or first not in START_CODONS
    three = (end - start + 1) % 3 != 0 or last not in STOP_CODONS
    return ("none", "5prime", "3prime", "both")[five + 2 * three]


def seq_hash(seq):
    return hashlib.sha1(seq.encode()).hexdigest()[:SHA1_CHARS]


def parse_gff_attrs(attr_str):
    out = {}
    for field in attr_str.strip().split(";"):
        field = field.strip()
        if field and "=" in field:
            k, v = field.split("=", 1)
            out[k] = v
    return out


def iter_orfipy(gff_path):
    """(orf_id, contig_key, start, end, strand, orf_len) from one genome's novel-ORF GFF."""
    with open(gff_path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                continue
            contig, _src, _type, start, end, _score, strand, _frame, attrs = cols
            a = parse_gff_attrs(attrs)
            yield (a.get("ID", ""), contig.rsplit("_", 1)[-1], int(start), int(end), strand,
                   a.get("orf_len"))


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--merged-table", required=True)
    ap.add_argument("--merged-summary", required=True)
    ap.add_argument("--orfipy-dir", required=True, type=Path)
    ap.add_argument("--genomes-dir", required=True, type=Path)
    ap.add_argument("--hygiene", required=True)
    ap.add_argument("--max-aa", type=int, default=100)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    gate = filterseq(yaml.safe_load(Path(args.hygiene).read_text()))

    with open(args.merged_summary) as f:
        next(f)
        genomes = [line.split("\t", 1)[0] for line in f]
    merged_rows = defaultdict(list)
    with open(args.merged_table) as f:
        next(f)
        for line in f:
            cols = line.rstrip("\n").split("\t")
            merged_rows[cols[0]].append(cols)

    per_genome = defaultdict(lambda: [0, 0, 0])   # gid -> [n_merged, n_orfipy, n_dropped]
    drop_reasons = defaultdict(int)
    length_hist = Counter()                        # (source, aa_length) -> n, for the summary
    partial_n = Counter()                          # (source, partial) -> n, for the summary
    n_len_mismatch = n_no_orfipy = 0
    p = args.prefix
    with open(f"{p}.smorfs.faa", "w") as faa, open(f"{p}.smorfs_metadata.tsv", "w") as meta, \
            open(f"{p}.dropped.tsv", "w") as drop:
        meta.write("genome_id\torf_id\tsource\tcontig\tstart\tend\tstrand"
                   "\taa_length\tpartial\ttools\tsupport_n\tseq_sha1\n")
        drop.write("genome_id\torf_id\tsource\treason\n")

        def dropped(gid, orf_id, source, reason):
            drop.write(f"{gid}\t{orf_id}\t{source}\t{reason}\n")
            drop_reasons[reason] += 1
            per_genome[gid][2] += 1

        for gid in genomes:
            contigs = load_contigs(args.genomes_dir / f"{gid}.fna")
            for cols in merged_rows.get(gid, []):
                _g, orf_id, contig, start, end, strand, aa_length, seq, tools, support_n, _l = cols
                reason = gate(seq)
                if reason:
                    dropped(gid, orf_id, "merged_catalog", reason)
                    continue
                part = partial_ends(contigs[contig], int(start), int(end), strand)
                faa.write(f">{orf_id}\n{seq}\n")
                meta.write(f"{gid}\t{orf_id}\tmerged_catalog\t{contig}\t{start}\t{end}\t{strand}"
                           f"\t{aa_length}\t{part}\t{tools}\t{support_n}\t{seq_hash(seq)}\n")
                per_genome[gid][0] += 1
                length_hist[("merged_catalog", int(aa_length))] += 1
                partial_n[("merged_catalog", part)] += 1

            gff = args.orfipy_dir / f"{gid}.gff"
            if not gff.exists():
                n_no_orfipy += 1
                continue
            for orf_id, contig_key, start, end, strand, orf_len in iter_orfipy(gff):
                if contig_key not in contigs:
                    print(f"WARN {gid}: contig {contig_key} not in FASTA, skipping {orf_id}",
                          file=sys.stderr)
                    continue
                protein = translate_nt(contigs[contig_key], start, end, strand)
                if len(protein) > args.max_aa:
                    print(f"WARN {gid}: {orf_id} translated to {len(protein)} aa > "
                          f"{args.max_aa}, skipping", file=sys.stderr)
                    continue
                if orf_len is not None and len(protein) != int(orf_len) // 3:
                    n_len_mismatch += 1
                reason = gate(protein)
                if reason:
                    dropped(gid, orf_id, "orfipy_novel", reason)
                    continue
                part = partial_ends(contigs[contig_key], start, end, strand)
                faa.write(f">{orf_id}\n{protein}\n")
                meta.write(f"{gid}\t{orf_id}\torfipy_novel\t{contig_key}\t{start}\t{end}"
                           f"\t{strand}\t{len(protein)}\t{part}\torfipy\t1\t{seq_hash(protein)}\n")
                per_genome[gid][1] += 1
                length_hist[("orfipy_novel", len(protein))] += 1
                partial_n[("orfipy_novel", part)] += 1

    with open(f"{p}.combine_map.tsv", "w") as m:
        m.write("genome_id\tn_merged\tn_orfipy\tn_dropped\tn_records\n")
        for gid in genomes:
            nm, no, nd = per_genome[gid]
            m.write(f"{gid}\t{nm}\t{no}\t{nd}\t{nm + no}\n")
    with open(f"{p}.combine_stats.tsv", "w") as s:
        s.write("kind\tkey\tvalue\tn\n")
        for (source, aa), n in sorted(length_hist.items()):
            s.write(f"aa_length\t{source}\t{aa}\t{n}\n")
        for reason, n in sorted(drop_reasons.items()):
            s.write(f"dropped\t{reason}\t\t{n}\n")
        for (source, part), n in sorted(partial_n.items()):
            s.write(f"partial\t{source}\t{part}\t{n}\n")

    reasons = " ".join(f"{k}={v}" for k, v in sorted(drop_reasons.items()))
    print(f"{p}: genomes={len(genomes)} merged_rows={sum(v[0] for v in per_genome.values())} "
          f"orfipy_rows={sum(v[1] for v in per_genome.values())} no_orfipy={n_no_orfipy} "
          f"len_mismatches={n_len_mismatch} dropped={sum(drop_reasons.values())} [{reasons}]")


if __name__ == "__main__":
    main()

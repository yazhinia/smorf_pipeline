#!/usr/bin/env python3
"""Drop orfipy ORFs that overlap an already-annotated gene, and emit GFF3.

orfipy is exhaustive: on a 2.7 Mbp genome it calls ~55k ORFs under the catalog's
parameters, against ~2.3k annotated genes. Only the ORFs in unannotated space are
interesting here, so the overlap filter runs in the worker and the unfiltered BED never
leaves node-local scratch.

Overlap is strand-agnostic and >= 1 bp: an ORF antisense to a called gene is not novel
coding space. See GENE_FEATURES for what counts as annotated; `region` is deliberately
excluded because it spans the whole contig and would drop every ORF.

    filter_novel_orfs.py --bed g.bed --merged-gff g.gff --genome-id g --out g.novel.gff
"""

import argparse
import os
import sys
from bisect import bisect_right
from functools import lru_cache

# normalize_seqname maps every tool's contig naming onto the bare `c001` token. The merged
# GFFs use `<genome>-c001` and orfipy's BED uses `<genome>_c001`; it handles both.
# merge_tools.py reaches into the same module for the same reason.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from merge_gffs import normalize_seqname as _normalize_seqname  # noqa: E402

# Called once per BED line -- ~55k times per genome against a few hundred distinct contig
# names, and the uncached regex was costing more than orfipy itself.
normalize_seqname = lru_cache(maxsize=None)(_normalize_seqname)

# Compared case-insensitively, so "sORF"/"sorf" and "crispr"/"CRISPR" both match.
# Only CDS/tRNA/rRNA/ncRNA/tmRNA occur in the merged GFFs today; the rest are bakta
# feature types that appear if --skip-crispr/--skip-gap are dropped, or come from other
# annotation sources.
# `region` spans the whole contig and `regulatory_region` is not a gene: both stay out.
GENE_FEATURES = {
    "cds",
    "gene",
    "pseudogene",
    "sorf",
    "trna",
    "rrna",
    "ncrna",
    "tmrna",
    "ncrna-region",
    "crispr",
    "gap",
    "is",
    "mite",
}


def load_gene_intervals(path):
    """contig -> (starts, ends) of disjoint 0-based half-open annotated intervals."""
    by_contig = {}
    n_features = 0
    with open(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2].lower() not in GENE_FEATURES:
                continue
            n_features += 1
            by_contig.setdefault(normalize_seqname(cols[0]), []).append(
                (int(cols[3]) - 1, int(cols[4]))
            )

    merged = {}
    for contig, ivs in by_contig.items():
        ivs.sort()
        starts, ends = [], []
        for a, b in ivs:
            if starts and a <= ends[-1]:
                ends[-1] = max(ends[-1], b)
            else:
                starts.append(a)
                ends.append(b)
        merged[contig] = (starts, ends)
    return merged, n_features


def overlaps(intervals, contig, a, b):
    """True if [a, b) touches any annotated interval on this contig."""
    entry = intervals.get(contig)
    if entry is None:
        return False
    starts, ends = entry
    # Intervals are disjoint and sorted, so only the last one starting before b can reach
    # back into [a, b).
    i = bisect_right(starts, b - 1) - 1
    return i >= 0 and ends[i] > a


def parse_orf_name(name):
    """orfipy BED name field -> (orf_len, frame, start_codon, stop_codon)."""
    orf_len = frame = start_codon = stop_codon = "."
    for field in name.split(";"):
        if field.startswith("ORF_len="):
            orf_len = field[8:]
        elif field.startswith("ORF_frame="):
            frame = field[10:]
        elif field.startswith("Start:"):
            start_codon = field[6:]
        elif field.startswith("Stop:"):
            stop_codon = field[5:]
    return orf_len, frame, start_codon, stop_codon


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--bed", required=True, help="orfipy BED (0-based half-open)")
    ap.add_argument("--merged-gff", required=True)
    ap.add_argument("--genome-id", required=True)
    ap.add_argument("--out", required=True, help="GFF3 of the surviving ORFs")
    ap.add_argument("--stats", help="TSV to append one summary line to")
    args = ap.parse_args()

    intervals, n_features = load_gene_intervals(args.merged_gff)

    n_raw = n_kept = 0
    with open(args.bed) as bed, open(args.out, "w") as out:
        out.write("##gff-version 3\n")
        for line in bed:
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 6:
                continue
            n_raw += 1
            contig = normalize_seqname(cols[0])
            a, b = int(cols[1]), int(cols[2])
            if overlaps(intervals, contig, a, b):
                continue
            n_kept += 1
            orf_len, frame, start_codon, stop_codon = parse_orf_name(cols[3])
            # Output seqid keeps the underscore form '<genome>_c001'; combine_chunk.py
            # splits genome and contig on the last '_'.
            out.write(
                f"{args.genome_id}_{contig}\torfipy\tCDS\t{a + 1}\t{b}\t.\t{cols[5]}\t0\t"
                f"ID={args.genome_id}_orfipy_{n_kept:06d};orf_len={orf_len};"
                f"start_codon={start_codon};stop_codon={stop_codon};frame={frame}\n"
            )

    if args.stats:
        # One line, one write: array elements append concurrently and single short
        # appends survive that, whereas a shared stdout log does not.
        with open(args.stats, "a") as fh:
            fh.write(f"{args.genome_id}\t{n_raw}\t{n_kept}\t{n_features}\n")

    print(f"{args.genome_id}: {n_raw} ORFs, {n_kept} kept, {n_features} gene features")


if __name__ == "__main__":
    main()

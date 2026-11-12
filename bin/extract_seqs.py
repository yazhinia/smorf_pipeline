#!/usr/bin/env python3
"""Unique protein and nucleotide sequences of one chunk, for AntiFam and Popcorn.

    extract_seqs.py --faa chunk.smorfs.faa --metadata chunk.smorfs_metadata.tsv \
        --genomes-dir genomes --prefix chunk_000001

The same smORF is found in many genomes, so each flag tool runs once per distinct sequence:
  <prefix>.aa_unique.faa   one record per seq_sha1 (AntiFam input)
  <prefix>.nt_unique.ffn   one record per distinct ORF nucleotide sequence (Popcorn input)
  <prefix>.nt_map.tsv      orf_id -> nt_id

Nucleotides include the stop codon, as extract_rep_nucleotide.py:93-98 wrote them for the
catalog run of Popcorn.
"""
import argparse
import hashlib
import sys
from collections import defaultdict
from pathlib import Path

from Bio.Seq import Seq

sys.path.insert(0, str(Path(__file__).resolve().parent))
from merge_chunk import load_contigs  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--faa", required=True)
    ap.add_argument("--metadata", required=True)
    ap.add_argument("--genomes-dir", required=True, type=Path)
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    # Proteins: the faa ids are orf_ids; metadata gives their seq_sha1.
    sha_of = {}
    by_genome = defaultdict(list)
    with open(args.metadata) as f:
        header = next(f).rstrip("\n").split("\t")
        ix = {c: i for i, c in enumerate(header)}
        for line in f:
            c = line.rstrip("\n").split("\t")
            sha_of[c[ix["orf_id"]]] = c[ix["seq_sha1"]]
            by_genome[c[ix["genome_id"]]].append(
                (c[ix["orf_id"]], c[ix["contig"]], int(c[ix["start"]]), int(c[ix["end"]]),
                 c[ix["strand"]]))

    n_aa = 0
    seen = set()
    with open(args.faa) as fin, open(f"{args.prefix}.aa_unique.faa", "w") as fout:
        orf_id = None
        for line in fin:
            if line.startswith(">"):
                orf_id = line[1:].split()[0]
                continue
            sha = sha_of[orf_id]
            if sha not in seen:
                seen.add(sha)
                fout.write(f">{sha}\n{line.strip()}\n")
                n_aa += 1

    nt_ids = {}
    with open(f"{args.prefix}.nt_unique.ffn", "w") as ffn, \
            open(f"{args.prefix}.nt_map.tsv", "w") as nmap:
        nmap.write("orf_id\tnt_id\n")
        for gid, rows in by_genome.items():
            contigs = load_contigs(args.genomes_dir / f"{gid}.fna")
            for orf_id, contig, start, end, strand in rows:
                nt = contigs[contig][start - 1:end]
                if strand == "-":
                    nt = str(Seq(nt).reverse_complement())
                if len(nt) % 3:
                    nt = nt[:len(nt) - len(nt) % 3]
                nid = hashlib.sha1(nt.encode()).hexdigest()[:20]
                if nid not in nt_ids:
                    nt_ids[nid] = None
                    ffn.write(f">{nid}\n{nt}\n")
                nmap.write(f"{orf_id}\t{nid}\n")

    print(f"{args.prefix}: orfs={len(sha_of)} unique_protein={n_aa} unique_nt={len(nt_ids)}")


if __name__ == "__main__":
    main()

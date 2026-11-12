#!/usr/bin/env python3
"""Write one chunk's genomes as normalized FASTAs with safe contig ids.

    normalize_chunk.py --plan plan.tsv --chunk-id chunk_000001 --split-mode file \
        [--split-regex RE] --prefix chunk_000001

Writes
  <prefix>.genomes/<genome_id>.fna   contigs renamed <genome_id>_c001.. (00_prepare_inputs.py:262)
  <prefix>.fna                       all genomes of the chunk concatenated (GMSC-Mapper input)
  <prefix>.contig_map.tsv            genome_id, sample_name, original_contig_id,
                                     normalized_contig_id, length
  <prefix>.genome_stats.tsv          one row per planned genome, including failed ones

Sequence text is copied as-is (no case change), exactly as the catalog's inputs were, and
empty records are dropped (00_prepare_inputs.py:231-233).
"""
import argparse
import gzip
import re
import shutil
import sys
from pathlib import Path

LINE_WIDTH = 80
STATS_HEADER = ("genome_id\tsample_name\tsource_fasta\tstatus\tdetail\tn_contigs\ttotal_bp"
                "\tn50\tgc_pct\tnon_acgt\n")


def open_text(path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path)


def iter_fasta(path):
    """(record_id, sequence) pairs; same parsing as 00_prepare_inputs.py:123-139."""
    header, chunks = None, []
    with open_text(path) as fh:
        for raw in fh:
            line = raw.strip()
            if not line:
                continue
            if line.startswith(">"):
                if header is not None:
                    yield header, "".join(chunks)
                header = line[1:].split()[0] if line[1:].strip() else ""
                chunks = []
            else:
                chunks.append(line)
    if header is not None:
        yield header, "".join(chunks)


def write_record(fh, name, seq):
    fh.write(f">{name}\n")
    for pos in range(0, len(seq), LINE_WIDTH):
        fh.write(seq[pos:pos + LINE_WIDTH] + "\n")


def n50(lengths):
    half, run = sum(lengths) / 2, 0
    for n in sorted(lengths, reverse=True):
        run += n
        if run >= half:
            return n
    return 0


class Genome:
    def __init__(self, row, out_dir):
        self.row = row
        self.gid = row["genome_id"]
        self.path = out_dir / f"{self.gid}.fna"
        self.fh = None
        self.lengths = []
        self.gc = self.non_acgt = 0
        self.contig_rows = []
        self.status, self.detail = "ok", ""

    def add(self, original_id, seq):
        if not seq:
            return
        if self.fh is None:
            self.fh = open(self.path, "w")
        name = f"{self.gid}_c{len(self.lengths) + 1:03d}"
        write_record(self.fh, name, seq)
        up = seq.upper()
        acgt = sum(up.count(b) for b in "ACGT")
        self.gc += up.count("G") + up.count("C")
        self.non_acgt += len(seq) - acgt
        self.lengths.append(len(seq))
        self.contig_rows.append((original_id, name, len(seq)))

    def close(self):
        if self.fh:
            self.fh.close()
        if self.status == "ok" and not self.lengths:
            self.status = "empty"
        if self.status != "ok":
            self.path.unlink(missing_ok=True)

    def stats_line(self):
        r = self.row
        total = sum(self.lengths)
        gc = f"{100 * self.gc / total:.2f}" if total else ""
        return (f"{self.gid}\t{r['sample_name']}\t{r['source_fasta']}\t{self.status}"
                f"\t{self.detail}\t{len(self.lengths)}\t{total}\t{n50(self.lengths)}\t{gc}"
                f"\t{self.non_acgt}\n")


def read_plan(path, chunk_id):
    rows = []
    with open(path) as fh:
        header = next(fh).rstrip("\n").split("\t")
        for line in fh:
            row = dict(zip(header, line.rstrip("\n").split("\t")))
            if row["chunk_id"] == chunk_id:
                rows.append(row)
    if not rows:
        sys.exit(f"no genomes for {chunk_id} in {path}")
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--plan", required=True)
    ap.add_argument("--chunk-id", required=True)
    ap.add_argument("--split-mode", choices=("file", "record", "regex"), default="file")
    ap.add_argument("--split-regex")
    ap.add_argument("--prefix", required=True)
    args = ap.parse_args()

    rows = read_plan(args.plan, args.chunk_id)
    out_dir = Path(f"{args.prefix}.genomes")
    out_dir.mkdir(parents=True, exist_ok=True)
    genomes = [Genome(r, out_dir) for r in rows]

    # Group by source file so a multi-genome FASTA is read once per chunk.
    by_file = {}
    for g in genomes:
        by_file.setdefault(g.row["source_fasta"], []).append(g)

    rx = re.compile(args.split_regex) if args.split_mode == "regex" else None
    for src, members in by_file.items():
        try:
            if args.split_mode == "file":
                for g in members:          # one genome per file; a file listed twice is read twice
                    for rid, seq in iter_fasta(src):
                        g.add(rid, seq)
            else:
                wanted = {g.row["selector"]: g for g in members}
                for rid, seq in iter_fasta(src):
                    name = rid if rx is None else rx.search(rid).group(1)
                    g = wanted.get(name)
                    if g is not None:
                        g.add(rid, seq)
        except (OSError, EOFError, UnicodeDecodeError, AttributeError) as exc:
            for g in members:
                g.status, g.detail = "unreadable", str(exc).replace("\t", " ")[:200]
        for g in members:
            g.close()

    with open(f"{args.prefix}.fna", "w") as chunk_fa, \
            open(f"{args.prefix}.contig_map.tsv", "w") as cmap, \
            open(f"{args.prefix}.genome_stats.tsv", "w") as stats:
        cmap.write("genome_id\tsample_name\toriginal_contig_id\tnormalized_contig_id\tlength\n")
        stats.write(STATS_HEADER)
        for g in genomes:
            stats.write(g.stats_line())
            if g.status != "ok":
                continue
            for orig, name, length in g.contig_rows:
                cmap.write(f"{g.gid}\t{g.row['sample_name']}\t{orig}\t{name}\t{length}\n")
            with open(g.path) as fh:
                shutil.copyfileobj(fh, chunk_fa)

    n_ok = sum(g.status == "ok" for g in genomes)
    print(f"{args.chunk_id}: genomes={len(genomes)} ok={n_ok} failed={len(genomes) - n_ok}")


if __name__ == "__main__":
    main()

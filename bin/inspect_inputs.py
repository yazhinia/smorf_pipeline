#!/usr/bin/env python3
"""Quick look at an input before a run: what it is and whether it is one genome or many.

    inspect_inputs.py <fasta|dir|samplesheet> [--max-files 50]

Reads headers and sequence lengths only, and stops early on big inputs. Prints a short
report plus a suggested split_mode, with the reason, and example genome/contig names (flagged
when long or with dots, spaces or '|', for the id_prefix question). Standard library only.
For a directory, only the first --max-files files are opened; the rest are just counted.
"""
import argparse
import gzip
import os
import re
import sys
from collections import Counter
from pathlib import Path

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")


def open_text(path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path)


def scan(path, max_records=200000):
    """-> list of (record_id, length), at most max_records."""
    recs, name, n = [], None, 0
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                if name is not None:
                    recs.append((name, n))
                    if len(recs) >= max_records:
                        return recs, True
                name, n = (line[1:].split() or [""])[0], 0
            else:
                n += len(line.strip())
    if name is not None:
        recs.append((name, n))
    return recs, False


def prefix_groups(ids):
    """Group record ids by a shared prefix such as 'genomeA_contig1' / 'genomeA|1'."""
    for sep in ("|", "__", "_", "."):
        keys = Counter(i.rsplit(sep, 1)[0] for i in ids if sep in i)
        if keys and sum(keys.values()) == len(ids) and 1 < len(keys) < len(ids):
            return sep, keys
    return None, None


def fasta_stem(path):
    lower = path.name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return path.name[: -len(suffix)]
    return path.stem


def name_note(label, names):
    """Example names, and how many are long or have characters other tools may cut or reject."""
    hard = [n for n in names if len(n) > 20 or not re.fullmatch(r"[A-Za-z0-9_-]+", n)]
    line = f"{label} (first 3): {', '.join(names[:3])}"
    if hard:
        line += (f"  [{len(hard)} of {len(names)} long or with dots/spaces/'|', "
                 f"e.g. {hard[0]}: short pipeline ids help]")
    print(line)


def describe_fasta(path):
    recs, truncated = scan(path)
    lengths = [n for _, n in recs]
    total = sum(lengths)
    print(f"file: {path}")
    print(f"records: {len(recs)}{'+ (stopped early)' if truncated else ''}  total_bp: {total:,}"
          f"  longest: {max(lengths, default=0):,}  shortest: {min(lengths, default=0):,}")
    ids = [r for r, _ in recs]
    print(f"first ids: {', '.join(ids[:5])}")
    name_note("genome name", [fasta_stem(path)])
    name_note("contig names", ids)
    big = sum(1 for n in lengths if n >= 500_000)
    if len(recs) == 1:
        print("suggest split_mode=file: one record, so one genome")
        return
    sep, groups = prefix_groups(ids)
    if groups:
        rx = r"^(.+)" + re.escape(sep) + r"[^" + re.escape(sep) + r"]+$"
        print(f"record ids share {len(groups)} prefixes before '{sep}' "
              f"(e.g. {', '.join(list(groups)[:3])})")
        print(f"suggest split_mode=regex --split_regex '{rx}' if each prefix is one genome; "
              f"otherwise split_mode=file")
    elif big >= 2 and big == len(recs):
        print(f"{big} records are each >= 500 kb, like whole genomes: "
              f"suggest split_mode=record if each record is a separate genome")
    else:
        print("looks like one assembly in several contigs: suggest split_mode=file "
              "(ask the user if unsure)")
    if total < 100_000:
        print("WARN: under 100 kb in total; is this really a genome?")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("input")
    ap.add_argument("--max-files", type=int, default=50)
    args = ap.parse_args()
    p = Path(args.input)
    if p.is_dir():
        files = sorted(Path(d) / f for d, _, fs in os.walk(p) for f in fs
                       if f.lower().endswith(FASTA_SUFFIXES))
        print(f"directory: {p}  FASTA files: {len(files)}  (one genome per file: split_mode=file)")
        if files:
            name_note("genome names", [fasta_stem(f) for f in files[:args.max_files]])
        sizes = []
        for f in files[:args.max_files]:
            recs, _ = scan(f, 100000)
            if not sizes:
                name_note("contig names in first file", [r for r, _ in recs])
            sizes.append((sum(n for _, n in recs), len(recs), f.name))
        if sizes:
            bp = sorted(s[0] for s in sizes)
            print(f"first {len(sizes)} files: genome size min/median/max = "
                  f"{bp[0]:,} / {bp[len(bp) // 2]:,} / {bp[-1]:,} bp; "
                  f"contigs median = {sorted(s[1] for s in sizes)[len(sizes) // 2]}")
        return
    if not p.is_file():
        sys.exit(f"not found: {p}")
    if p.suffix.lower() in (".tsv", ".csv"):
        delim = "," if p.suffix.lower() == ".csv" else "\t"
        lines = p.read_text().splitlines()
        header = lines[0].split(delim)
        print(f"samplesheet: {p}  rows: {len(lines) - 1}  columns: {', '.join(header)}")
        missing = {"sample", "fasta"} - set(header)
        if missing:
            print(f"ERROR: missing column(s) {', '.join(sorted(missing))}")
        else:
            col = header.index("sample")
            name_note("genome names", [ln.split(delim)[col] for ln in lines[1:] if ln.strip()])
        return
    describe_fasta(p)


if __name__ == "__main__":
    main()

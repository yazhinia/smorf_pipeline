#!/usr/bin/env python3
"""Decide which genomes a run has, give each a safe id, and assign chunks.

    plan_inputs.py --input <fasta|dir|samplesheet> --split-mode file|record|regex \
        [--split-regex RE] --id-prefix gU --chunk-size 500 --out plan.tsv

Input forms:
  - one FASTA file             split_mode decides whether it is one genome or many
  - a directory                every FASTA under it (recursive), one genome per file
  - a samplesheet (.tsv/.csv)  columns sample, fasta and optionally genome_id

FASTA discovery and suffixes follow 00_prepare_inputs.py (:15-22, :50-86); the HRGM2
metadata join is dropped. In file mode no FASTA is opened, so planning 100k genomes is fast.
In record / regex mode the headers are read once here to know the genomes.

plan.tsv columns: genome_id, sample_name, source_fasta, selector, chunk_id
  selector is '*' (whole file) or the sample name whose records belong to this genome.
"""
import argparse
import csv
import gzip
import os
import re
import sys
from collections import Counter
from pathlib import Path

FASTA_SUFFIXES = (".fa", ".fna", ".fasta", ".fa.gz", ".fna.gz", ".fasta.gz")
# No '_', '-' or '.': merge_gffs.py splits file ids on '-' and the contig/ORF ids are
# '<genome>_c001' / '<genome>_ORF00001', which are split on '_' downstream.
GENOME_ID_RE = re.compile(r"^[A-Za-z][A-Za-z0-9]{2,15}$")


def fasta_stem(path):
    lower = path.name.lower()
    for suffix in FASTA_SUFFIXES:
        if lower.endswith(suffix):
            return path.name[: -len(suffix)]
    return path.stem


def is_fasta(path):
    return path.name.lower().endswith(FASTA_SUFFIXES)


def open_text(path):
    return gzip.open(path, "rt") if str(path).lower().endswith(".gz") else open(path)


def record_ids(path):
    with open_text(path) as fh:
        for line in fh:
            if line.startswith(">"):
                yield line[1:].split()[0] if line[1:].strip() else ""


def read_samplesheet(path):
    delim = "," if path.suffix.lower() == ".csv" else "\t"
    rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter=delim)
        missing = {"sample", "fasta"} - set(reader.fieldnames or [])
        if missing:
            sys.exit(f"samplesheet {path} is missing column(s): {', '.join(sorted(missing))}")
        for n, row in enumerate(reader, start=2):
            fasta = Path(row["fasta"].strip())
            if not fasta.is_absolute():
                fasta = (path.parent / fasta).resolve()
            if not fasta.is_file():
                sys.exit(f"samplesheet line {n}: FASTA not found: {fasta}")
            rows.append({"sample": row["sample"].strip(), "fasta": fasta,
                         "genome_id": (row.get("genome_id") or "").strip()})
    if not rows:
        sys.exit(f"samplesheet {path} has no rows")
    return rows


def discover(input_path):
    """-> list of {sample, fasta, genome_id} for file-level inputs, sorted by path."""
    p = Path(input_path)
    if p.is_dir():
        found = sorted(Path(d) / f for d, _, files in os.walk(p.absolute())
                       for f in files if f.lower().endswith(FASTA_SUFFIXES))
        if not found:
            sys.exit(f"no FASTA files ({', '.join(FASTA_SUFFIXES)}) under {p}")
        return [{"sample": fasta_stem(f), "fasta": f, "genome_id": ""} for f in found]
    if not p.is_file():
        sys.exit(f"input not found: {p}")
    if p.suffix.lower() in (".tsv", ".csv"):
        return read_samplesheet(p)
    if not is_fasta(p):
        sys.exit(f"unsupported input {p}: expected FASTA ({', '.join(FASTA_SUFFIXES)}), "
                 f"a directory, or a .tsv/.csv samplesheet")
    return [{"sample": fasta_stem(p), "fasta": p.absolute(), "genome_id": ""}]


def expand(entries, split_mode, split_regex):
    """file mode: one genome per entry. record/regex mode: genomes from FASTA record ids."""
    if split_mode == "file":
        return [dict(e, selector="*") for e in entries]
    rx = re.compile(split_regex) if split_mode == "regex" else None
    if rx is not None and rx.groups < 1:
        sys.exit(f"--split-regex needs one capture group: {split_regex}")
    genomes = []
    for e in entries:
        seen = {}
        for rid in record_ids(e["fasta"]):
            if rx is None:
                name = rid
            else:
                m = rx.search(rid)
                if not m:
                    sys.exit(f"record '{rid}' in {e['fasta']} does not match --split-regex")
                name = m.group(1)
            if name not in seen:
                seen[name] = None
        for name in seen:
            genomes.append({"sample": name, "fasta": e["fasta"], "genome_id": "",
                            "selector": name})
    return genomes


def assign_ids(genomes, prefix):
    explicit = [g["genome_id"] for g in genomes if g["genome_id"]]
    if explicit and len(explicit) != len(genomes):
        sys.exit("samplesheet: give genome_id for every row or for none")
    if not explicit and not GENOME_ID_RE.match(f"{prefix}000001"):
        sys.exit(f"--id-prefix '{prefix}' must be letters/digits and start with a letter")
    for n, g in enumerate(genomes, start=1):
        if not explicit:
            g["genome_id"] = f"{prefix}{n:06d}"
        elif not GENOME_ID_RE.match(g["genome_id"]):
            sys.exit(f"genome_id '{g['genome_id']}' must match {GENOME_ID_RE.pattern}")
    dup = [i for i, n in Counter(g["genome_id"] for g in genomes).items() if n > 1]
    if dup:
        sys.exit(f"duplicate genome_id(s): {sorted(dup)[:5]}")
    names = [g["sample"] for g in genomes]
    if len(set(names)) != len(names):
        print("WARN: duplicate sample names; genome_id stays unique", file=sys.stderr)


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--input", required=True)
    ap.add_argument("--split-mode", choices=("file", "record", "regex"), default="file")
    ap.add_argument("--split-regex")
    ap.add_argument("--id-prefix", default="gU")
    ap.add_argument("--chunk-size", type=int, default=500)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()
    if args.split_mode == "regex" and not args.split_regex:
        sys.exit("--split-mode regex needs --split-regex")
    if args.chunk_size < 1:
        sys.exit("--chunk-size must be >= 1")

    genomes = expand(discover(args.input), args.split_mode, args.split_regex)
    if not genomes:
        sys.exit("no genomes found")
    assign_ids(genomes, args.id_prefix)

    with open(args.out, "w") as out:
        out.write("genome_id\tsample_name\tsource_fasta\tselector\tchunk_id\n")
        for n, g in enumerate(genomes):
            chunk_id = f"chunk_{n // args.chunk_size + 1:06d}"
            out.write(f"{g['genome_id']}\t{g['sample']}\t{g['fasta']}\t{g['selector']}"
                      f"\t{chunk_id}\n")

    n_chunks = (len(genomes) - 1) // args.chunk_size + 1
    print(f"genomes={len(genomes)} chunks={n_chunks} split_mode={args.split_mode}")


if __name__ == "__main__":
    main()

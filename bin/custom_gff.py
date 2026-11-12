#!/usr/bin/env python3
"""Clean one genome's GFF from an extra caller for the merge, and record its status.

    custom_gff.py --tool NAME --genome gU000001.fna --raw out.gff --rc 0 \
        --out gU000001.gff --status status.tsv

Keeps CDS rows on the genome's own contigs (<gid>_cNNN), writes the tool name into column 2
(the merge reads the caller from there), and gives every row an ID. A tool that exits non-zero
or writes no GFF is a soft failure: one status row, no GFF, the task goes on.
Standard library only, so it runs under any tool environment.
"""
import argparse
import sys
from pathlib import Path


def contig_ids(fna):
    with open(fna) as f:
        return {line[1:].split()[0] for line in f if line.startswith(">")}


def clean(raw, tool, contigs):
    """-> (kept GFF lines, n CDS on unknown contigs, first unknown seqid)."""
    kept, n_unknown, first_unknown = [], 0, ""
    with open(raw) as f:
        for line in f:
            if line.startswith("##FASTA"):
                break
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9 or cols[2].lower() != "cds" or cols[6] not in "+-":
                continue
            seqid = cols[0].split()[0]
            if seqid not in contigs:
                n_unknown += 1
                first_unknown = first_unknown or seqid
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                continue
            if start > end:
                continue
            attrs = cols[8].strip().rstrip(";")
            if not any(kv.startswith("ID=") for kv in attrs.split(";")):
                attrs = f"ID={tool}_{len(kept) + 1}" + (f";{attrs}" if attrs and attrs != "." else "")
            kept.append("\t".join([seqid, tool, "CDS", str(start), str(end),
                                   cols[5], cols[6], cols[7], attrs]))
    return kept, n_unknown, first_unknown


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--tool", required=True)
    ap.add_argument("--genome", required=True, type=Path)
    ap.add_argument("--raw", required=True, type=Path)
    ap.add_argument("--rc", type=int, default=0)
    ap.add_argument("--out", required=True, type=Path)
    ap.add_argument("--status", required=True, type=Path)
    args = ap.parse_args()

    gid = args.genome.name.rsplit(".", 1)[0]
    status, reason, n = "failed", "", 0
    if args.rc == 124:
        reason = "timeout"
    elif args.rc:
        reason = f"rc_{args.rc}"
    elif not args.raw.exists():
        reason = "no_gff"
    else:
        kept, n_unknown, first_unknown = clean(args.raw, args.tool, contig_ids(args.genome))
        if n_unknown and not kept:
            reason = f"unknown_contigs:{first_unknown}"
        else:
            n = len(kept)
            status = "ok" if kept else "empty_ok"
            reason = f"unknown_contig_rows={n_unknown}" if n_unknown else ("" if kept else "no_cds")
            with open(args.out, "w") as fh:
                fh.write("##gff-version 3\n")
                fh.writelines(line + "\n" for line in kept)

    new = not args.status.exists()
    with open(args.status, "a") as fh:
        if new:
            fh.write("genome_id\ttool\tstatus\treason\tn_records\n")
        fh.write(f"{gid}\t{args.tool}\t{status}\t{reason}\t{n}\n")
    if status == "failed":
        print(f"{gid} {args.tool} FAILED {reason}", file=sys.stderr)


if __name__ == "__main__":
    main()

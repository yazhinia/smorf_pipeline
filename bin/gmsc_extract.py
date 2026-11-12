#!/usr/bin/env python3
"""Turn one GMSC-Mapper chunk output into per-genome GFFs of high-quality smORFs.

    gmsc_extract.py --gmsc-out <gmsc-mapper -o dir> --contig-map chunk.contig_map.tsv \
        --out-dir gmsc/

Port of process_gmscmapper_chunks.sh:12-17:
  1. smORFs whose quality.out.smorfs.tsv row says 'high quality' and that are present in
     mapped.smorfs.faa
  2. their headers in predicted.filterd.smorf.faa, '>ID # contig_orf # start # end # strand # attrs',
     become GFF lines 'contig_orf gmscmapper CDS start end . +/- . Name=ID;attrs'
  3. one <genome_id>.gff per genome that has at least one such smORF

The shell version matched ids with grep -F (substring); this matches whole ids, which is
identical while ids are fixed-width (smORF_NNNNNN) and stays correct beyond that.
Genomes come from the contig map rather than from splitting seqids on '_c'.
"""
import argparse
import sys
from pathlib import Path


def fasta_ids(path):
    with open(path) as fh:
        return {line[1:].split()[0] for line in fh if line.startswith(">")}


def high_quality_ids(path):
    hq = set()
    with open(path) as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 2 and cols[1] == "high quality":
                hq.add(cols[0])
    return hq


def contig_to_genome(path):
    m = {}
    with open(path) as fh:
        next(fh)
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            m[cols[3]] = cols[0]
    return m


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gmsc-out", required=True, type=Path)
    ap.add_argument("--contig-map", required=True)
    ap.add_argument("--out-dir", required=True, type=Path)
    args = ap.parse_args()

    d = args.gmsc_out
    keep = high_quality_ids(d / "quality.out.smorfs.tsv") & fasta_ids(d / "mapped.smorfs.faa")
    genome_of = contig_to_genome(args.contig_map)
    args.out_dir.mkdir(parents=True, exist_ok=True)

    handles, counts = {}, {}
    with open(d / "predicted.filterd.smorf.faa") as fh:
        for line in fh:
            if not line.startswith(">"):
                continue
            a = line[1:].rstrip("\n").split(" # ")
            if a[0] not in keep:
                continue
            a += [""] * (6 - len(a))
            seqid = a[1]                                  # <genome>_c001_<orf index>
            contig = seqid.rsplit("_", 1)[0]
            gid = genome_of.get(contig)
            if gid is None:
                sys.exit(f"GMSC contig {contig} (from {seqid}) is not in the contig map")
            if gid not in handles:
                handles[gid] = open(args.out_dir / f"{gid}.gff", "w")
                handles[gid].write("##gff-version 3\n")
                counts[gid] = 0
            strand = "+" if a[4] == "1" else "-"
            handles[gid].write(f"{seqid}\tgmscmapper\tCDS\t{a[2]}\t{a[3]}\t.\t{strand}\t.\t"
                               f"Name={a[0]};{a[5]}\n")
            counts[gid] += 1
    for h in handles.values():
        h.close()

    with open(args.out_dir / "status.tsv", "w") as st:
        st.write("genome_id\ttool\tstatus\treason\tn_records\n")
        for gid in dict.fromkeys(genome_of.values()):
            st.write(f"{gid}\tgmscmapper\tok\t\t{counts.get(gid, 0)}\n")
    print(f"high-quality smORFs kept: {sum(counts.values())} across {len(counts)} genomes")


if __name__ == "__main__":
    main()

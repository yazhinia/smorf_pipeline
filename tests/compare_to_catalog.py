#!/usr/bin/env python3
"""Check that the pipeline reproduces the microprotein catalog for a catalog genome.

    # merge logic only, from the catalog's own per-tool GFFs (no tools run)
    compare_to_catalog.py fast-merge --genome gI053058 --workdir DIR

    # a finished pipeline run (-profile test, publish_intermediates = true)
    compare_to_catalog.py run --genome gI053058 --results RUN/results

Catalog files are located through its batch maps, so any catalog genome id works.
Set SMORF_CATALOG_DIR to the catalog's orf_prediction_run folder (and SMORF_CATALOG_NFS_DIR
to its NFS copy, if separate). Prints one PASS/FAIL line per check and exits 1 if any fails.
"""
import argparse
import csv
import glob
import os
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bin"))

LUSTRE = Path(os.environ.get("SMORF_CATALOG_DIR", "orf_prediction_run"))
NFS = Path(os.environ.get("SMORF_CATALOG_NFS_DIR", str(LUSTRE)))
TOOLS = LUSTRE / "tool_outputs"
NFS_TOOLS = NFS / "tool_outputs"

failures = []


def check(name, ok, detail=""):
    print(f"{'PASS' if ok else 'FAIL'}  {name}{'  ' + detail if detail else ''}")
    if not ok:
        failures.append(name)


def map_lookup(path, gid, col):
    with open(path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["genome_id"] == gid:
                return r[col]
    sys.exit(f"{gid} not in {path}")


def rows_for(path, gid, key="genome_id"):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return [r for r in reader if r[key] == gid or r[key].startswith(f"{gid}_")]


def gff_lines(path, prefix, ncols=9):
    out = []
    with open(path) as fh:
        for line in fh:
            if line.startswith(prefix):
                out.append(tuple(line.rstrip("\n").split("\t")[:ncols]))
    return sorted(out)


def catalog_sources(gid):
    bakta_file = map_lookup(TOOLS / "bakta_combined/bakta_gff3_batch_map.tsv", gid, "combined_file") \
        if (TOOLS / "bakta_combined/bakta_gff3_batch_map.tsv").exists() else None
    smorf_file = map_lookup(NFS_TOOLS / "smorfinder_combined/smorfinder_gff_batch_map.tsv",
                            gid, "combined_gff")
    with open(NFS / "indexes/gmsc_chunk_index.tsv") as fh:
        chunk = next(r["chunk_id"] for r in csv.DictReader(fh, delimiter="\t")
                     if r["genome_id"] == gid)
    batch = map_lookup(TOOLS / "merged_orf_table_batch_map.tsv", gid, "batch_id")
    return {
        "bakta": TOOLS / "bakta_combined" / bakta_file if bakta_file else TOOLS / f"bakta/{gid}.gff3",
        "smorf": NFS_TOOLS / "smorfinder_combined" / smorf_file,
        "gmsc": NFS_TOOLS / f"gmscmapper/{chunk}/{gid}.gff",
        "batch": batch,
        "fasta": NFS / f"inputs_normalized/{gid}.fna",
    }


def fast_merge(args):
    from merge_chunk import genome_feature_lines, load_contigs, translate
    from merge_gffs import merge_gff_files

    src = catalog_sources(args.genome)
    wd = Path(args.workdir)
    wd.mkdir(parents=True, exist_ok=True)
    gid = args.genome
    bakta = wd / f"{gid}.gff3"
    with open(src["bakta"]) as fin, open(bakta, "w") as fout:
        fout.write("##gff-version 3\n")
        for line in fin:
            if line.startswith(f"{gid}_c") or line.startswith(f"##sequence-region {gid}_c"):
                fout.write(line)
    smorf = wd / "smorf" / f"{gid}.gff"
    smorf.parent.mkdir(exist_ok=True)
    genome_feature_lines(src["smorf"], smorf, gid)
    inputs = [str(bakta)] + [str(p) for p in (src["gmsc"], smorf) if p.exists()]
    kept = merge_gff_files(inputs, str(wd / f"{gid}.merged.gff"))
    contigs = load_contigs(src["fasta"])
    ours = []
    for orf_id, contig, start, end, strand, tools in kept:
        protein = translate(contigs[contig], start, end, strand)
        if len(protein) <= 100:
            ours.append((gid, orf_id, contig, str(start), str(end), strand, str(len(protein)),
                         protein, ",".join(tools), str(len(tools)), ""))
    ref = [tuple(r.values()) for r in rows_for(TOOLS / f"merged_orf_table/{src['batch']}.tsv", gid)]
    check("merged_orf_table rows identical", sorted(ours) == sorted(ref),
          f"ours={len(ours)} catalog={len(ref)}")


def run(args):
    gid = args.genome
    res = Path(args.results)
    src = catalog_sources(gid)
    chunk_dir = next(Path(p) for p in glob.glob(str(res / "chunks/chunk_*"))
                     if rows_for(next(iter(glob.glob(f"{p}/*.genome_stats.tsv"))), gid))
    cid = chunk_dir.name
    inter = res / "intermediates"

    # 1. callers (need publish_intermediates = true)
    b = next(iter(glob.glob(str(inter / f"{cid}.bakta_b*/{gid}.gff3"))), None)
    if b:
        ours = gff_lines(b, f"{gid}_c")
        ref = gff_lines(src["bakta"], f"{gid}_c")
        check("bakta GFF (all columns)", ours == ref, f"ours={len(ours)} catalog={len(ref)}")
        if ours != ref:
            o8 = sorted(x[:8] for x in ours)
            r8 = sorted(x[:8] for x in ref)
            check("bakta GFF (columns 1-8)", o8 == r8, f"diff={len(set(o8) ^ set(r8))}")
    s = next(iter(glob.glob(str(inter / f"{cid}.smorf_b*/{gid}.gff"))), None)
    if s:
        ours = sorted(x[:8] for x in gff_lines(s, f"{gid}_c"))
        ref = sorted(x[:8] for x in gff_lines(src["smorf"], f"{gid}_c"))
        check("smORFinder GFF (columns 1-8)", ours == ref, f"ours={len(ours)} catalog={len(ref)}")
    g = inter / f"{cid}.gmsc/{gid}.gff"
    if g.exists():
        ours = sorted(x[:8] for x in gff_lines(g, f"{gid}_c"))
        ref = sorted(x[:8] for x in gff_lines(src["gmsc"], f"{gid}_c"))
        check("GMSC-Mapper GFF (columns 1-8)", ours == ref, f"ours={len(ours)} catalog={len(ref)}")

    # 2. merge (later stages are skipped when a run stopped at an earlier checkpoint)
    batch = src["batch"]
    t = chunk_dir / f"{cid}.merged_orf_table.tsv"
    if t.exists():
        ours = sorted(tuple(r.values()) for r in rows_for(t, gid))
        ref = sorted(tuple(r.values()) for r in rows_for(TOOLS / f"merged_orf_table/{batch}.tsv", gid))
        check("merged_orf_table rows identical", ours == ref, f"ours={len(ours)} catalog={len(ref)}")
    if not (chunk_dir / f"{cid}.merged_summary.tsv").exists():
        print("SKIP  merge, combine and flags: not run yet")
        return
    ours = rows_for(chunk_dir / f"{cid}.merged_summary.tsv", gid)
    ref = rows_for(TOOLS / f"merged_summary/{batch}.tsv", gid)
    check("merged_summary counts identical",     # catalog files may carry later columns
          [tuple(r.values())[:7] for r in ours] == [tuple(r.values())[:7] for r in ref],
          f"ours={[list(r.values())[1:] for r in ours]} catalog={[list(r.values())[1:7] for r in ref]}")

    # 3. combine (combined_batch numbering follows the merged batch)
    if not (chunk_dir / f"{cid}.smorfs.tsv").exists():
        print("SKIP  combine and flags: not run yet")
        return
    meta_cols = ["genome_id", "orf_id", "source", "contig", "start", "end", "strand",
                 "aa_length", "tools", "support_n", "seq_sha1"]
    smorfs = rows_for(chunk_dir / f"{cid}.smorfs.tsv", gid)
    ours = sorted(tuple(r[c] for c in meta_cols) for r in smorfs)
    ref_meta = TOOLS / f"combined_microproteins/combined_{batch}_metadata.tsv"
    ref = sorted(tuple(r[c] for c in meta_cols) for r in rows_for(ref_meta, gid))
    check("combined smORFs identical (11 metadata columns)", ours == ref,
          f"ours={len(ours)} catalog={len(ref)} only_ours={len(set(ours) - set(ref))} "
          f"only_catalog={len(set(ref) - set(ours))}")

    # 4. flags, where the catalog has them (it flagged cluster representatives only)
    ours_by_id = {r["orf_id"]: r for r in smorfs}
    pc_ref = {}
    pc_path = TOOLS / "microprotein_catalog_clusters/popcorn/popcorn_output.tsv"
    with open(pc_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["Sequence ID"].startswith(f"{gid}_"):
                pc_ref[r["Sequence ID"]] = r
    shared = [i for i in pc_ref if i in ours_by_id]
    bad = [i for i in shared
           if abs(float(ours_by_id[i]["popcorn_probability"] or -9)
                  - float(pc_ref[i]["Coding Probability"])) > 1e-4]
    check("Popcorn probability on catalog reps", bool(shared) and not bad,
          f"compared={len(shared)} differ={len(bad)}")
    af_path = LUSTRE / "classifier_dataset/evidence/antifam_flags.tsv"
    af_ref = {}
    with open(af_path) as fh:
        for r in csv.DictReader(fh, delimiter="\t"):
            if r["rep_orf_id"].startswith(f"{gid}_"):
                af_ref[r["rep_orf_id"]] = r["antifam_hit"]
    shared = [i for i in af_ref if i in ours_by_id]
    bad = [i for i in shared if ours_by_id[i]["antifam_hit"] != af_ref[i]]
    check("AntiFam hit flag on catalog reps", bool(shared) and not bad,
          f"compared={len(shared)} differ={len(bad)}")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    sub = ap.add_subparsers(dest="mode", required=True)
    f = sub.add_parser("fast-merge")
    f.add_argument("--genome", required=True)
    f.add_argument("--workdir", required=True)
    r = sub.add_parser("run")
    r.add_argument("--genome", required=True)
    r.add_argument("--results", required=True)
    args = ap.parse_args()
    if not LUSTRE.is_dir():
        sys.exit(f"catalog folder {LUSTRE} not found: set SMORF_CATALOG_DIR")
    fast_merge(args) if args.mode == "fast-merge" else run(args)
    print(f"\n{'ALL PASS' if not failures else f'{len(failures)} FAILED'}")
    sys.exit(1 if failures else 0)


if __name__ == "__main__":
    main()

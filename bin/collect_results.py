#!/usr/bin/env python3
"""Final run-level tables from the per-chunk outputs (staged in the working directory).

    collect_results.py --concat-max-genomes 1000

Always writes
  per_genome_summary.tsv  one row per planned genome: input stats, merge and combine counts
  genome_map.tsv          genome_id -> sample name and source FASTA
  contig_map.tsv          normalized contig id -> original contig name
Up to --concat-max-genomes genomes, also writes the smORFs as one set of files
  smorfs.tsv, smorfs.faa, smorfs.gff3 (on the user's original contig names), dropped.tsv
Above that they stay per chunk (results/chunks/<chunk>/), since a 100k-genome run is
~1.7 billion rows.
"""
import argparse
import glob


def read_tsv(path):
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def concat(pattern, out_path):
    """Concatenate same-header TSVs, keeping the first header only."""
    n = 0
    with open(out_path, "w") as out:
        for i, p in enumerate(sorted(glob.glob(pattern))):
            with open(p) as f:
                header = next(f)
                if i == 0:
                    out.write(header)
                for line in f:
                    out.write(line)
                    n += 1
    return n


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--concat-max-genomes", type=int, default=1000)
    args = ap.parse_args()

    stats = [r for p in sorted(glob.glob("*.genome_stats.tsv")) for r in read_tsv(p)]
    merged = {r["genome_id"]: r for p in glob.glob("*.merged_summary.tsv") for r in read_tsv(p)}
    skipped = {r["genome_id"]: r["reason"]
               for p in glob.glob("*.merge_skipped.tsv") for r in read_tsv(p)}
    combined = {r["genome_id"]: r for p in glob.glob("*.combine_map.tsv") for r in read_tsv(p)}

    with open("genome_map.tsv", "w") as gm:
        gm.write("genome_id\tsample_name\tsource_fasta\tstatus\n")
        for r in stats:
            gm.write(f"{r['genome_id']}\t{r['sample_name']}\t{r['source_fasta']}\t{r['status']}\n")
    concat("*.contig_map.tsv", "contig_map.tsv")

    cols = ["genome_id", "sample_name", "input_status", "n_contigs", "total_bp", "gc_pct",
            "merge_status", "n_cds", "n_microprotein_merged", "n_smorfs_merged_catalog",
            "n_smorfs_orfipy_novel", "n_dropped_hygiene"]
    with open("per_genome_summary.tsv", "w") as out:
        out.write("\t".join(cols) + "\n")
        for r in stats:
            g = r["genome_id"]
            m, c = merged.get(g, {}), combined.get(g, {})
            if r["status"] != "ok":
                mstatus = "not_run"
            elif g in merged:
                mstatus = "merged"
            else:
                mstatus = skipped.get(g, "not_run")
            vals = [g, r["sample_name"], r["status"], r["n_contigs"], r["total_bp"],
                    r["gc_pct"], mstatus, m.get("n_cds", ""), m.get("n_microprotein", ""),
                    c.get("n_merged", ""), c.get("n_orfipy", ""), c.get("n_dropped", "")]
            out.write("\t".join(vals) + "\n")

    n_ok = sum(r["status"] == "ok" for r in stats)
    if n_ok > args.concat_max_genomes:
        print(f"{n_ok} genomes > {args.concat_max_genomes}: smORF tables stay per chunk")
        return

    n = concat("*.smorfs.tsv", "smorfs.tsv")
    concat("*.dropped.tsv", "dropped.tsv")
    with open("smorfs.faa", "w") as faa, open("smorfs.gff3", "w") as gff:
        gff.write("##gff-version 3\n")
        for r in read_tsv("smorfs.tsv"):
            faa.write(f">{r['orf_id']} genome={r['genome_id']} sample={r['sample_name']} "
                      f"source={r['source']}\n{r['sequence']}\n")
            attrs = ";".join([
                f"ID={r['orf_id']}", f"genome_id={r['genome_id']}",
                f"smorf_source={r['source']}", f"tools={r['tools']}",
                f"support_n={r['support_n']}", f"aa_length={r['aa_length']}",
                f"partial={r['partial']}",
                f"antifam_hit={r['antifam_hit']}",
                f"popcorn_probability={r['popcorn_probability']}",
            ])
            gff.write(f"{r['original_contig']}\tsmorf_pipeline\tCDS\t{r['start']}\t{r['end']}"
                      f"\t.\t{r['strand']}\t0\t{attrs}\n")
    print(f"smorfs={n} genomes_ok={n_ok}")


if __name__ == "__main__":
    main()

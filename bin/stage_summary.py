#!/usr/bin/env python3
"""Checkpoint summary for one pipeline stage, read by the user (and the smorf-predict skill).

    stage_summary.py --stage prepare|call|merge|combine|flags --out-dir qc [inputs in cwd]

Reads only the small per-chunk files each stage writes (stats, statuses, counts), never the
full ORF tables, so it stays quick on 100k-genome runs. Writes
  qc/<stage>_summary.tsv      metric, value, level (OK/WARN/FAIL), note
  qc/<stage>_summary.md       <= 20 lines for the checkpoint conversation
  qc/<stage>_per_genome.tsv   one row per genome
  qc/<stage>_*.png            combine and flags only
"""
import argparse
import glob
import statistics
from collections import Counter, defaultdict
from pathlib import Path

BAKTA_SEC_PER_GENOME = 1238   # mean wall time at 4 threads in the catalog run (logs/bakta)


def read_tsv(path):
    with open(path) as f:
        header = next(f).rstrip("\n").split("\t")
        for line in f:
            if line.strip():
                yield dict(zip(header, line.rstrip("\n").split("\t")))


def rows(pattern):
    for p in sorted(glob.glob(pattern)):
        yield from read_tsv(p)


def pct(a, b):
    return f"{100 * a / b:.1f}" if b else "0.0"


class Summary:
    def __init__(self, stage):
        self.stage = stage
        self.metrics = []
        self.per_genome = []
        self.per_genome_cols = []
        self.figures = []

    def add(self, metric, value, level="OK", note=""):
        self.metrics.append((metric, str(value), level, note))

    def write(self, out_dir):
        out = Path(out_dir)
        out.mkdir(parents=True, exist_ok=True)
        with open(out / f"{self.stage}_summary.tsv", "w") as f:
            f.write("metric\tvalue\tlevel\tnote\n")
            for m in self.metrics:
                f.write("\t".join(m) + "\n")
        if self.per_genome_cols:
            with open(out / f"{self.stage}_per_genome.tsv", "w") as f:
                f.write("\t".join(self.per_genome_cols) + "\n")
                for r in self.per_genome:
                    f.write("\t".join(str(r.get(c, "")) for c in self.per_genome_cols) + "\n")
        levels = Counter(m[2] for m in self.metrics)
        status = "FAIL" if levels["FAIL"] else ("WARN" if levels["WARN"] else "OK")
        lines = [f"## Stage `{self.stage}`: {status}", "",
                 "| metric | value | level |", "|---|---|---|"]
        lines += [f"| {m} | {v} | {lv} |" for m, v, lv, _ in self.metrics[:12]]
        notes = [f"- **{m}**: {n}" for m, _, lv, n in self.metrics if lv != "OK" and n]
        lines += [""] + notes[:4] if notes else []
        if self.figures:
            lines += ["", "Figures: " + ", ".join(self.figures)]
        (out / f"{self.stage}_summary.md").write_text("\n".join(lines[:20]) + "\n")


def setup_matplotlib():
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams.update({
        "font.family": "sans-serif",
        "font.sans-serif": ["Arial", "Helvetica", "Calibri", "DejaVu Sans"],
        "font.size": 7, "axes.labelsize": 7, "axes.titlesize": 8, "legend.fontsize": 6,
        "xtick.labelsize": 6, "ytick.labelsize": 6, "axes.linewidth": 0.75,
        "lines.linewidth": 1.25, "figure.autolayout": True,
    })
    return plt


SINGLE_COLUMN = 85 / 25.4
GOLDEN_RATIO = 1.618
SOURCE_COLOURS = {"merged_catalog": "#1f6fb4", "orfipy_novel": "#d9822b"}
BUILTIN_TOOLS = ("bakta", "smorfinder", "gmscmapper")


def prepare(s):
    stats = list(rows("*.genome_stats.tsv"))
    ok = [r for r in stats if r["status"] == "ok"]
    failed = [r for r in stats if r["status"] != "ok"]
    s.add("genomes_planned", len(stats))
    s.add("genomes_ok", len(ok), "FAIL" if not ok else "OK",
          "no readable genome: check the input path and FASTA format" if not ok else "")
    s.add("genomes_failed", len(failed), "WARN" if failed else "OK",
          f"empty or unreadable: {', '.join(r['sample_name'] for r in failed[:5])}"
          if failed else "")
    if ok:
        sizes = [int(r["total_bp"]) for r in ok]
        contigs = [int(r["n_contigs"]) for r in ok]
        total_bp = sum(sizes)
        non_acgt = sum(int(r["non_acgt"]) for r in ok)
        many = [r for r in ok if int(r["n_contigs"]) > 999]
        s.add("contigs_total", sum(contigs))
        s.add("genome_size_bp_min_median_max",
              f"{min(sizes)} / {int(statistics.median(sizes))} / {max(sizes)}")
        s.add("gc_pct_median", f"{statistics.median(float(r['gc_pct']) for r in ok):.1f}")
        s.add("non_acgt_pct", pct(non_acgt, total_bp),
              "WARN" if non_acgt > 0.01 * total_bp else "OK",
              "over 1% of bases are not A/C/G/T (N runs or IUPAC codes)")
        s.add("genomes_over_999_contigs", len(many), "WARN" if many else "OK",
              "contig c1000 sorts before c101; only ORF numbering order changes")
        n_chunks = len(set(glob.glob("*.genome_stats.tsv")))
        s.add("chunks", n_chunks)
        s.add("est_bakta_core_hours", f"{len(ok) * BAKTA_SEC_PER_GENOME * 4 / 3600:.0f}",
              note="Bakta is the slowest step (~21 min per genome on 4 cores)")
        s.add("gmsc_200gb_jobs", n_chunks,
              note="one 16-core, 200 GB GMSC-Mapper job per chunk, even for one genome")
    s.per_genome_cols = ["genome_id", "sample_name", "status", "detail", "n_contigs",
                         "total_bp", "n50", "gc_pct", "non_acgt", "source_fasta"]
    s.per_genome = stats


def call(s):
    by_tool = defaultdict(dict)
    for p in sorted(glob.glob("*/status.tsv")):
        for r in read_tsv(p):
            by_tool[r["tool"]][r["genome_id"]] = r
    versions = {}                     # extra callers with version_cmd write <dir>/version.txt
    for p in sorted(glob.glob("*/version.txt")):
        tool, _, v = Path(p).read_text().rstrip("\n").partition("\t")
        versions[tool] = v
    extra = sorted(t for t in by_tool if t not in BUILTIN_TOOLS)
    tools = list(BUILTIN_TOOLS) + extra
    genomes = sorted({g for t in by_tool.values() for g in t})
    s.add("tools", ", ".join(tools),
          note="gene callers run; extra callers merge after the built-in three" if extra else "")
    for tool in tools:
        res = by_tool.get(tool, {})
        n_ok = sum(r["status"] in ("ok", "empty_ok") for r in res.values())
        n_fail = len(res) - n_ok
        reasons = Counter(r["reason"] for r in res.values() if r["status"] == "failed")
        level = "OK"
        if n_fail:
            level = "FAIL" if tool not in ("smorfinder", "gmscmapper") and n_ok == 0 else "WARN"
        note = ", ".join(f"{k}={v}" for k, v in reasons.most_common(3))
        if tool == "bakta" and n_fail:
            note += "; genomes without Bakta are skipped by the merge"
        elif tool in extra and n_fail:
            note += "; these genomes get no calls from this tool"
        s.add(f"{tool}_ok", f"{n_ok}/{len(res)}", level, note)
        if tool in versions:
            s.add(f"{tool}_version", versions[tool])
        recs = [int(r["n_records"]) for r in res.values()
                if r["status"] != "failed" and r["n_records"].isdigit()]
        if recs:
            s.add(f"{tool}_records_per_genome_median", int(statistics.median(recs)))
    s.per_genome_cols = ["genome_id"] + [f"{t}_{c}" for t in tools
                                         for c in ("status", "reason", "n_records")]
    for g in genomes:
        row = {"genome_id": g}
        for t in tools:
            r = by_tool.get(t, {}).get(g, {})
            for c in ("status", "reason", "n_records"):
                row[f"{t}_{c}"] = r.get(c, "missing")
        s.per_genome.append(row)


def merge(s):
    summ = list(rows("*.merged_summary.tsv"))
    skipped = list(rows("*.merge_skipped.tsv"))
    combos = Counter()
    for r in rows("*.merge_combos.tsv"):
        combos[r["tools"]] += int(r["n_microprotein"])
    n_micro = sum(int(r["n_microprotein"]) for r in summ)
    n_stop = sum(int(r["n_internal_stop"]) for r in summ)
    s.add("genomes_merged", len(summ), "FAIL" if not summ else "OK")
    s.add("genomes_skipped", len(skipped), "WARN" if skipped else "OK",
          "no Bakta GFF for these genomes" if skipped else "")
    s.add("cds_total", sum(int(r["n_cds"]) for r in summ))
    s.add("microproteins_total", n_micro, note="merged calls <= max_aa")
    if summ:
        s.add("microproteins_per_genome_median",
              int(statistics.median(int(r["n_microprotein"]) for r in summ)))
    s.add("internal_stop", n_stop, "WARN" if n_stop else "OK",
          "translated with an internal stop; removed at the combine stage")
    for tools, n in combos.most_common(7):
        s.add(f"support:{tools}", f"{n} ({pct(n, n_micro)}%)")
    s.per_genome_cols = list(summ[0]) if summ else ["genome_id", "n_cds", "n_microprotein"]
    s.per_genome = summ


def combine(s, out_dir):
    cmap = list(rows("*.combine_map.tsv"))
    hist = defaultdict(Counter)
    dropped = Counter()
    partial = defaultdict(Counter)     # source -> {5prime, 3prime, both, none: n}
    for r in rows("*.combine_stats.tsv"):
        if r["kind"] == "aa_length":
            hist[r["key"]][int(r["value"])] += int(r["n"])
        elif r["kind"] == "dropped":
            dropped[r["key"]] += int(r["n"])
        elif r["kind"] == "partial":
            partial[r["key"]][r["value"]] += int(r["n"])
    n_m = sum(int(r["n_merged"]) for r in cmap)
    n_o = sum(int(r["n_orfipy"]) for r in cmap)
    n_d = sum(dropped.values())
    no_orfipy = sum(1 for r in cmap if r["n_orfipy"] == "0")
    s.add("smorfs_merged_catalog", n_m, note="called by Bakta, smORFinder and/or GMSC-Mapper")
    s.add("smorfs_orfipy_novel", n_o, note="ORFs in unannotated space; unlabelled discovery pool")
    s.add("dropped_hygiene", n_d, "WARN" if n_d > 0.01 * max(n_m + n_o, 1) else "OK",
          ", ".join(f"{k}={v}" for k, v in dropped.most_common()))
    for src, n_src in (("merged_catalog", n_m), ("orfipy_novel", n_o)):
        c = partial.get(src, Counter())
        n_p = n_src - c["none"]
        s.add(f"partial:{src}", n_p, "WARN" if n_p > 0.01 * max(n_src, 1) else "OK",
              ", ".join(f"{k}={c[k]}" for k in ("5prime", "3prime", "both") if c[k]))
    s.add("genomes_without_orfipy_rows", no_orfipy, "WARN" if no_orfipy else "OK",
          "orfipy found no ORF outside annotated genes" if no_orfipy else "")
    for src, h in sorted(hist.items()):
        half, run = sum(h.values()) / 2, 0
        for aa in sorted(h):
            run += h[aa]
            if run >= half:
                s.add(f"aa_length_median:{src}", aa)
                break
    s.per_genome_cols = ["genome_id", "n_merged", "n_orfipy", "n_dropped", "n_records"]
    s.per_genome = cmap
    if hist:
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(SINGLE_COLUMN, SINGLE_COLUMN / GOLDEN_RATIO))
        for src, h in sorted(hist.items()):
            xs = list(range(min(h), max(h) + 1))     # lengths with no smORF plot as 0
            total = sum(h.values())
            ax.plot(xs, [100 * h[x] / total for x in xs], label=f"{src} (n={total:,})",
                    color=SOURCE_COLOURS.get(src))
        ax.set_xlabel("Protein length (aa)")
        ax.set_ylabel("smORFs (% of source)")
        ax.legend(frameon=False)
        ax.margins(x=0.02)
        name = "combine_aa_length.png"
        fig.savefig(Path(out_dir) / name, dpi=600, bbox_inches="tight", pad_inches=0.04)
        s.figures.append(f"{name} (x: protein length in amino acids; y: share of each "
                         f"source's smORFs at that length)")


def flags(s, out_dir):
    counts = defaultdict(Counter)
    bins = defaultdict(Counter)
    for r in rows("*.flags_stats.tsv"):
        if r["kind"] == "count":
            counts[r["key"]][r["value"]] += int(r["n"])
        elif r["kind"] == "popcorn_bin":
            bins[r["key"]][float(r["value"])] += int(r["n"])
    for src, c in sorted(counts.items()):
        n = c["smorfs"]
        af = 100 * c["antifam_hit"] / n if n else 0
        miss = 100 * c["popcorn_missing"] / n if n else 0
        s.add(f"antifam_hit_pct:{src}", f"{af:.2f}", "WARN" if af > 25 else "OK",
              "over 25% flagged spurious: check input and HMM file")
        s.add(f"popcorn_coding_pct:{src}", pct(c["popcorn_coding"], n - c["popcorn_missing"]))
        s.add(f"popcorn_missing_pct:{src}", f"{miss:.2f}", "WARN" if miss > 5 else "OK",
              "Popcorn skipped these sequences; left empty, not called non-coding")
    if bins:
        plt = setup_matplotlib()
        fig, ax = plt.subplots(figsize=(SINGLE_COLUMN, SINGLE_COLUMN / GOLDEN_RATIO))
        for src, h in sorted(bins.items()):
            total = sum(h.values())
            xs = [i / 20 for i in range(20)]           # empty bins plot as 0
            ax.step([x + 0.025 for x in xs], [100 * h[x] / total for x in xs], where="mid",
                    label=f"{src} (n={total:,})", color=SOURCE_COLOURS.get(src))
        ax.set_xlabel("Popcorn coding probability")
        ax.set_ylabel("smORFs (% of source)")
        ax.set_xlim(0, 1)
        ax.legend(frameon=False)
        name = "flags_popcorn.png"
        fig.savefig(Path(out_dir) / name, dpi=600, bbox_inches="tight", pad_inches=0.04)
        s.figures.append(f"{name} (x: Popcorn coding probability in 0.05 bins; y: share of "
                         f"each source's smORFs in that bin)")


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--stage", required=True,
                    choices=("prepare", "call", "merge", "combine", "flags"))
    ap.add_argument("--out-dir", default="qc")
    args = ap.parse_args()
    s = Summary(args.stage)
    Path(args.out_dir).mkdir(parents=True, exist_ok=True)
    if args.stage == "prepare":
        prepare(s)
    elif args.stage == "call":
        call(s)
    elif args.stage == "merge":
        merge(s)
    elif args.stage == "combine":
        combine(s, args.out_dir)
    else:
        flags(s, args.out_dir)
    s.write(args.out_dir)
    print(Path(args.out_dir, f"{args.stage}_summary.md").read_text())


if __name__ == "__main__":
    main()

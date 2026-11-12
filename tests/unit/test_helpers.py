"""Unit tests for bin/ helpers on tiny synthetic inputs.

    python -m unittest discover -s tests/unit      (on a compute node, test_tools python)
"""
import gzip
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parents[2] / "bin"
sys.path.insert(0, str(BIN))

from hygiene import filterseq  # noqa: E402

# 300 nt contig: an ATG...TAA ORF of 7 codons at 11-31 on '+', plus filler.
ORF = "ATGAAACCCGGGTTTAAATAA"          # M K P G F K * (21 nt, stop included)
CONTIG = "GGGGGGGGGG" + ORF + "C" * 269


def run(script, *args, cwd):
    res = subprocess.run([sys.executable, str(BIN / script), *map(str, args)], cwd=cwd,
                         capture_output=True, text=True)
    if res.returncode:
        raise AssertionError(f"{script} failed:\n{res.stdout}\n{res.stderr}")
    return res.stdout


def tsv(path):
    lines = Path(path).read_text().splitlines()
    header = lines[0].split("\t")
    return [dict(zip(header, line.split("\t"))) for line in lines[1:]]


class TestHygiene(unittest.TestCase):
    def test_reasons(self):
        gate = filterseq({"hygiene": {"ambig_drop": {
            "drop_any_x": True, "drop_internal_stop": True, "nonstandard_residues": "BZJUO"}}})
        self.assertIsNone(gate("MKP"))
        self.assertIsNone(gate("MKP*"))
        self.assertEqual(gate("MXP"), "ambiguous_x")
        self.assertEqual(gate("MK*P"), "internal_stop")
        self.assertEqual(gate("MKB"), "nonstandard_residue")
        self.assertEqual(gate(""), "empty")


class TestPrepare(unittest.TestCase):
    def test_file_and_record_modes(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "in").mkdir()
            (d / "in/a.fna").write_text(">x1 desc\nACGT\nACGT\n>x2\nGGCC\n")
            with gzip.open(d / "in/b.fa.gz", "wt") as fh:
                fh.write(">y1\nAAAA\n>empty\n\n")
            run("plan_inputs.py", "--input", d / "in", "--split-mode", "file",
                "--chunk-size", "1", "--out", "plan.tsv", cwd=d)
            plan = tsv(d / "plan.tsv")
            self.assertEqual([r["genome_id"] for r in plan], ["gU000001", "gU000002"])
            self.assertEqual([r["chunk_id"] for r in plan], ["chunk_000001", "chunk_000002"])
            self.assertEqual(plan[0]["sample_name"], "a")

            run("normalize_chunk.py", "--plan", "plan.tsv", "--chunk-id", "chunk_000001",
                "--prefix", "chunk_000001", cwd=d)
            fa = (d / "chunk_000001.genomes/gU000001.fna").read_text()
            self.assertEqual(fa, ">gU000001_c001\nACGTACGT\n>gU000001_c002\nGGCC\n")
            cmap = tsv(d / "chunk_000001.contig_map.tsv")
            self.assertEqual(cmap[0]["original_contig_id"], "x1")
            stats = tsv(d / "chunk_000001.genome_stats.tsv")
            self.assertEqual((stats[0]["status"], stats[0]["n_contigs"], stats[0]["total_bp"]),
                             ("ok", "2", "12"))

            # one multi-genome FASTA, grouped by record-name prefix
            (d / "multi.fna").write_text(">gA|1\nAC\n>gB|1\nGG\n>gA|2\nTT\n")
            run("plan_inputs.py", "--input", d / "multi.fna", "--split-mode", "regex",
                "--split-regex", r"^(.+)\|", "--out", "plan2.tsv", cwd=d)
            plan2 = tsv(d / "plan2.tsv")
            self.assertEqual([r["sample_name"] for r in plan2], ["gA", "gB"])
            run("normalize_chunk.py", "--plan", "plan2.tsv", "--chunk-id", "chunk_000001",
                "--split-mode", "regex", "--split-regex", r"^(.+)\|", "--prefix", "m", cwd=d)
            self.assertEqual((d / "m.genomes/gU000001.fna").read_text(),
                             ">gU000001_c001\nAC\n>gU000001_c002\nTT\n")

    def test_bad_genome_id(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "g.fna").write_text(">c\nACGT\n")
            (d / "s.tsv").write_text(f"sample\tfasta\tgenome_id\ns1\t{d / 'g.fna'}\tbad_id\n")
            res = subprocess.run([sys.executable, str(BIN / "plan_inputs.py"), "--input",
                                  str(d / "s.tsv"), "--out", str(d / "p.tsv")],
                                 capture_output=True, text=True)
            self.assertNotEqual(res.returncode, 0)


class TestGmscExtract(unittest.TestCase):
    def test_high_quality_only(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            raw = d / "raw"
            raw.mkdir()
            (raw / "quality.out.smorfs.tsv").write_text(
                "qseqid\tquality\nsmORF_000001\thigh quality\nsmORF_000002\tlow quality\n"
                "smORF_000003\thigh quality\n")
            (raw / "mapped.smorfs.faa").write_text(
                ">smORF_000001\nMKP*\n>smORF_000002\nMAA*\n")
            (raw / "predicted.filterd.smorf.faa").write_text(
                ">smORF_000001 # gU000001_c002_7 # 11 # 31 # -1 # ID=2_7;partial=00\nMKP*\n"
                ">smORF_000002 # gU000001_c001_1 # 1 # 30 # 1 # ID=1_1;partial=00\nMAA*\n"
                ">smORF_000003 # gU000001_c001_2 # 1 # 30 # 1 # ID=1_2;partial=00\nMCC*\n")
            (d / "cmap.tsv").write_text(
                "genome_id\tsample_name\toriginal_contig_id\tnormalized_contig_id\tlength\n"
                "gU000001\ts\tx1\tgU000001_c001\t300\ngU000001\ts\tx2\tgU000001_c002\t300\n")
            run("gmsc_extract.py", "--gmsc-out", raw, "--contig-map", "cmap.tsv",
                "--out-dir", "gmsc", cwd=d)
            lines = (d / "gmsc/gU000001.gff").read_text().splitlines()
            self.assertEqual(lines, [
                "##gff-version 3",
                "gU000001_c002_7\tgmscmapper\tCDS\t11\t31\t.\t-\t.\tName=smORF_000001;ID=2_7;partial=00",
            ])


class TestMergeCombineFlags(unittest.TestCase):
    def test_pipeline_on_one_genome(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            gid = "gU000001"
            (d / "genomes").mkdir()
            (d / f"genomes/{gid}.fna").write_text(f">{gid}_c001\n{CONTIG}\n")
            (d / "stats.tsv").write_text(
                "genome_id\tsample_name\tsource_fasta\tstatus\tdetail\tn_contigs\ttotal_bp\tn50"
                f"\tgc_pct\tnon_acgt\n{gid}\ts1\t/x\tok\t\t1\t300\t300\t50\t0\n")
            (d / "bakta").mkdir()
            (d / f"bakta/{gid}.gff3").write_text(
                f"##gff-version 3\n##sequence-region {gid}_c001 1 300\n"
                f"{gid}_c001\tPyrodigal\tCDS\t11\t31\t.\t+\t0\tID=a\n")
            (d / "smorf").mkdir()
            # same stop codon (end 31), later start: same gene, tool tag carries over
            (d / f"smorf/{gid}.gff").write_text(
                f"{gid}_c001\tProdigal_v2.6.3\tCDS\t14\t31\t5\t+\t0\tID=b\n")
            (d / "gmsc").mkdir()
            run("merge_chunk.py", "--genome-stats", "stats.tsv", "--genomes-dir", "genomes",
                "--bakta-dirs", "bakta", "--smorf-dirs", "smorf", "--gmsc-dir", "gmsc",
                "--prefix", "c1", cwd=d)
            table = tsv(d / "c1.merged_orf_table.tsv")
            self.assertEqual(len(table), 1)
            self.assertEqual((table[0]["sequence"], table[0]["tools"], table[0]["support_n"]),
                             ("MKPGFK", "bakta,smorfinder", "2"))
            self.assertEqual(table[0]["orf_id"], f"{gid}_ORF00001")

            (d / "orfipy").mkdir()
            (d / f"orfipy/{gid}.gff").write_text(
                f"##gff-version 3\n{gid}_c001\torfipy\tCDS\t200\t220\t.\t+\t0\t"
                f"ID={gid}_orfipy_000001;orf_len=18\n")
            (d / "hyg.yaml").write_text(
                "hygiene:\n  ambig_drop:\n    drop_any_x: true\n    drop_internal_stop: true\n"
                "    nonstandard_residues: BZJUO\n")
            run("combine_chunk.py", "--merged-table", "c1.merged_orf_table.tsv",
                "--merged-summary", "c1.merged_summary.tsv", "--orfipy-dir", "orfipy",
                "--genomes-dir", "genomes", "--hygiene", "hyg.yaml", "--prefix", "c1", cwd=d)
            meta = tsv(d / "c1.smorfs_metadata.tsv")
            self.assertEqual([r["source"] for r in meta], ["merged_catalog", "orfipy_novel"])
            self.assertEqual(meta[1]["aa_length"], "7")      # CCC.. -> PPPPPPP, no stop codon
            self.assertEqual([r["partial"] for r in meta], ["none", "both"])

            run("extract_seqs.py", "--faa", "c1.smorfs.faa", "--metadata",
                "c1.smorfs_metadata.tsv", "--genomes-dir", "genomes", "--prefix", "c1", cwd=d)
            self.assertEqual((d / "c1.aa_unique.faa").read_text().count(">"), 2)
            ffn = (d / "c1.nt_unique.ffn").read_text().splitlines()
            self.assertIn(ORF, ffn)

            nt_ids = {r["orf_id"]: r["nt_id"] for r in tsv(d / "c1.nt_map.tsv")}
            (d / "af.tbl").write_text(
                f"{meta[0]['seq_sha1']} - AntiFam_1 ANF00001.1 1e-10 40.0 0.0\n# [ok]\n")
            (d / "pc.tsv").write_text(
                "Sequence ID\tCoding Probability\tPrediction\n"
                f"{nt_ids[meta[0]['orf_id']]}\t0.91\tCODING\n")
            (d / "cmap.tsv").write_text(
                "genome_id\tsample_name\toriginal_contig_id\tnormalized_contig_id\tlength\n"
                f"{gid}\ts1\tchrom\t{gid}_c001\t300\n")
            run("join_flags.py", "--metadata", "c1.smorfs_metadata.tsv", "--faa", "c1.smorfs.faa",
                "--nt-map", "c1.nt_map.tsv", "--antifam", "af.tbl", "--popcorn", "pc.tsv",
                "--contig-map", "cmap.tsv", "--prefix", "c1", cwd=d)
            out = tsv(d / "c1.smorfs.tsv")
            self.assertEqual((out[0]["antifam_hit"], out[0]["popcorn_prediction"],
                              out[0]["original_contig"], out[0]["sample_name"]),
                             ("1", "CODING", "chrom", "s1"))
            self.assertEqual((out[1]["antifam_hit"], out[1]["popcorn_probability"]), ("0", ""))

            for stage in ("merge", "combine", "flags"):
                run("stage_summary.py", "--stage", stage, "--out-dir", "qc", cwd=d)
            self.assertTrue((d / "qc/combine_aa_length.png").exists())
            self.assertIn("microproteins_total", (d / "qc/merge_summary.tsv").read_text())
            self.assertIn("partial:orfipy_novel\t1", (d / "qc/combine_summary.tsv").read_text())

    def test_partial_ends(self):
        from Bio.Seq import Seq
        from combine_chunk import partial_ends
        rc = str(Seq(ORF).reverse_complement())
        contig = "GGG" + ORF + "GGG" + rc + "GGG"     # ORF at 4-24 (+) and 28-48 (-)
        cases = {(4, 24, "+"): "none",
                 (28, 48, "-"): "none",
                 (7, 24, "+"): "5prime",     # AAA.., no start codon
                 (4, 21, "+"): "3prime",     # ..AAA, no stop codon
                 (31, 48, "-"): "3prime",
                 (7, 21, "+"): "both"}
        for (start, end, strand), want in cases.items():
            self.assertEqual(partial_ends(contig, start, end, strand), want, (start, end, strand))
        # a full ORF in the contig's first codon (+) or last codon (-) may start off the contig
        self.assertEqual(partial_ends(ORF + "GGG", 1, 21, "+"), "5prime")
        self.assertEqual(partial_ends("GGG" + rc, 4, 24, "-"), "5prime")

    def test_truncated_antifam_fails(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "af.tbl").write_text("x - A B 1 2 0\n")
            res = subprocess.run([sys.executable, str(BIN / "join_flags.py"), "--metadata", "m",
                                  "--faa", "f", "--nt-map", "n", "--antifam", "af.tbl",
                                  "--popcorn", "p", "--contig-map", "c", "--prefix", "x"],
                                 cwd=d, capture_output=True, text=True)
            self.assertNotEqual(res.returncode, 0)
            self.assertIn("[ok]", res.stderr)


class TestLogDecision(unittest.TestCase):
    def test_one_log_readable_first(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            run("log_decision.py", "--run-dir", "r", "--stage", "setup", "--question",
                "Bacterial?", "--choice", "Bacterial", "--comment", "E. coli\tisolate", cwd=d)
            rows = tsv(d / "r/run_log/decisions.tsv")
            self.assertEqual(list(rows[0])[:6], ["timestamp", "stage", "question", "choice",
                                                 "comment", "params_diff"])
            self.assertEqual(rows[0]["comment"], "E. coli isolate")
            self.assertEqual(sorted(p.name for p in (d / "r/run_log").iterdir()),
                             ["decisions.tsv"])

    def test_old_column_order_kept(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            old = ["timestamp", "user", "run_dir", "stage", "question", "options", "choice",
                   "comment", "params_diff", "params_md5", "pipeline_git"]
            (d / "r/run_log").mkdir(parents=True)
            (d / "r/run_log/decisions.tsv").write_text("\t".join(old) + "\n")
            run("log_decision.py", "--run-dir", "r", "--stage", "call", "--question", "Go?",
                "--choice", "Stop here", cwd=d)
            row = tsv(d / "r/run_log/decisions.tsv")[0]
            self.assertEqual((row["stage"], row["choice"]), ("call", "Stop here"))


class TestExtraCallers(unittest.TestCase):
    gid = "gU000001"

    def genome(self, d):
        (d / "genomes").mkdir(exist_ok=True)
        fna = d / f"genomes/{self.gid}.fna"
        fna.write_text(f">{self.gid}_c001\n{CONTIG}\n")
        return fna

    def test_custom_gff(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fna = self.genome(d)
            (d / "raw.gff").write_text(
                "##gff-version 3\n"
                f"{self.gid}_c001\tProdigal_v2.6.3\tgene\t11\t31\t.\t+\t.\tID=g1\n"
                f"{self.gid}_c001\tProdigal_v2.6.3\tCDS\t11\t31\t.\t+\t0\tID=c1;partial=00\n"
                f"{self.gid}_c001 desc\tx\tcds\t40\t60\t.\t-\t0\t.\n"
                f"other_contig\tx\tCDS\t1\t30\t.\t+\t0\tID=c3\n"
                "##FASTA\n>x\nACGT\n")
            run("custom_gff.py", "--tool", "mytool", "--genome", fna, "--raw", "raw.gff",
                "--out", "out.gff", "--status", "status.tsv", cwd=d)
            rows = [line.split("\t") for line in (d / "out.gff").read_text().splitlines()[1:]]
            self.assertEqual([(r[1], r[3], r[8]) for r in rows],
                             [("mytool", "11", "ID=c1;partial=00"), ("mytool", "40", "ID=mytool_2")])
            st = tsv(d / "status.tsv")[0]
            self.assertEqual((st["status"], st["reason"], st["n_records"]),
                             ("ok", "unknown_contig_rows=1", "2"))
            run("custom_gff.py", "--tool", "mytool", "--genome", fna, "--raw", "none.gff",
                "--rc", "124", "--out", "x.gff", "--status", "status.tsv", cwd=d)
            self.assertEqual(tsv(d / "status.tsv")[1]["reason"], "timeout")
            self.assertFalse((d / "x.gff").exists())

    def test_run_caller_then_merge(self):
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            fna = self.genome(d)
            # same stop as Bakta's call below (end 31), later start: one gene, two callers
            (d / "cmd.txt").write_text(
                "c=$(head -1 {fasta} | cut -c2-)\n"
                "printf '%s\\tanything\\tCDS\\t14\\t31\\t.\\t+\\t0\\tID=q\\n' \"$c\" > {gff}\n")
            res = subprocess.run(["bash", str(BIN / "run_caller.sh"), "--name", "mytool",
                                  "--env", "none", "--command-file", "cmd.txt",
                                  "--out-dir", "c1.x_mytool_b0", "--python", sys.executable,
                                  "--version-cmd", "echo mytool 1.2", str(fna)],
                                 cwd=d, capture_output=True, text=True)
            self.assertEqual(res.returncode, 0, res.stderr)
            self.assertEqual(tsv(d / "c1.x_mytool_b0/status.tsv")[0]["status"], "ok")
            self.assertEqual((d / "c1.x_mytool_b0/version.txt").read_text(), "mytool\tmytool 1.2\n")

            (d / "stats.tsv").write_text(
                "genome_id\tsample_name\tsource_fasta\tstatus\tdetail\tn_contigs\ttotal_bp\tn50"
                f"\tgc_pct\tnon_acgt\n{self.gid}\ts1\t/x\tok\t\t1\t300\t300\t50\t0\n")
            (d / "bakta").mkdir()
            (d / f"bakta/{self.gid}.gff3").write_text(
                f"{self.gid}_c001\tPyrodigal\tCDS\t11\t31\t.\t+\t0\tID=a\n")
            (d / "bakta/status.tsv").write_text(
                f"genome_id\ttool\tstatus\treason\tn_records\n{self.gid}\tbakta\tok\t\t1\n")
            (d / "gmsc").mkdir()
            run("merge_chunk.py", "--genome-stats", "stats.tsv", "--genomes-dir", "genomes",
                "--bakta-dirs", "bakta", "--gmsc-dir", "gmsc", "--prefix", "c1",
                "--extra-tools", "mytool", "--extra-dirs", "c1.x_mytool_b0", cwd=d)
            table = tsv(d / "c1.merged_orf_table.tsv")
            self.assertEqual((table[0]["tools"], table[0]["support_n"], table[0]["start"]),
                             ("bakta,mytool", "2", "11"))
            summ = tsv(d / "c1.merged_summary.tsv")[0]
            self.assertEqual((summ["n_2tool"], summ["n_4tool"]), ("1", "0"))

            run("stage_summary.py", "--stage", "call", "--out-dir", "qc", cwd=d)
            text = (d / "qc/call_summary.tsv").read_text()
            self.assertIn("tools\tbakta, smorfinder, gmscmapper, mytool", text)
            self.assertIn("mytool_ok\t1/1\tOK", text)
            self.assertIn("mytool_version\tmytool 1.2", text)
            run("stage_summary.py", "--stage", "merge", "--out-dir", "qc", cwd=d)
            self.assertIn("n_4tool", (d / "qc/merge_per_genome.tsv").read_text())

    def test_check_callers(self):
        from list_tools import check_callers
        ok = {"name": "my_tool", "env": "none", "command": "x -i {fasta} -o {gff}"}
        self.assertEqual(check_callers([ok]), [])
        bad = [
            ({**ok, "name": "pyrodigal"}, "taken by a built-in"),
            ({**ok, "name": "My"}, "lowercase"),
            ({**ok, "command": "x -i {fasta}"}, "{gff}"),
            ({**ok, "command": "x -i '{fasta}' -o {gff}"}, "do not quote"),
            ({**ok, "env": "docker:x"}, "env must be"),
            ({**ok, "batch": 0}, "whole number"),
            ({**ok, "time": "2 hours"}, "'90m'"),
        ]
        for entry, msg in bad:
            errs = check_callers([entry])
            self.assertTrue(any(msg in e for e in errs), (entry, errs))
        self.assertTrue(any("used twice" in e for e in check_callers([ok, ok])))

    def test_tool_rows_all_stages(self):
        from list_tools import tool_rows
        with tempfile.TemporaryDirectory() as d:
            d = Path(d)
            (d / "env/bin").mkdir(parents=True)
            (d / "env/conda-meta").mkdir()
            (d / "env/conda-meta/orfipy-0.0.4-py312_4.json").write_text("{}")
            (d / "antifam").mkdir()
            (d / "antifam/version").write_text("AntiFam release: 8.0\nAntiFam families: 278\n")
            rows = tool_rows({"orfipy_bin": str(d / "env/bin/orfipy"),
                              "antifam_hmm": str(d / "antifam/AntiFam.hmm"),
                              "popcorn_dir": str(d), "orf_min": 15, "orf_max": 300,
                              "extra_callers": [{"name": "mytool", "env": "none",
                                                 "command": "x {fasta} {gff}"}]})
        by_tool = {r["tool"]: r for r in rows}
        self.assertEqual([r["stage"] for r in rows],
                         ["call"] * 4 + ["combine", "flags", "flags"])
        self.assertEqual(by_tool["orfipy"]["software"], "orfipy 0.0.4")
        self.assertEqual(by_tool["antifam"]["database"], "AntiFam release 8.0, 278 families")
        self.assertEqual(by_tool["popcorn"]["software"], "Popcorn (local copy)")


if __name__ == "__main__":
    unittest.main()

#!/usr/bin/env python3
"""Every prediction tool a run uses, by stage, and a check of its extra_callers.

    list_tools.py <run_dir> [--tsv FILE]        table by category: smORF prediction (gene callers,
                                                + extra), intergenic smORF, coding potential,
                                                spurious tag
    list_tools.py <run_dir> --shell NAME        one extra caller as shell variables (test_caller.sh)

Reads the pipeline's defaults (nextflow.config `params {}`, then conf/site.config if present)
and the run's params.yaml, the same
values Nextflow will use. Versions come from the installs (module names, conda metadata,
release files, git), without running the tools. Exit 2 if an extra_callers entry is invalid,
with the reason. Needs PyYAML (the pipeline's python_bin has it).
"""
import argparse
import re
import shlex
import subprocess
import sys
from pathlib import Path

import yaml

PIPE = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(PIPE / "bin"))
from merge_gffs import SOURCE_TO_TOOL  # noqa: E402

BUILTIN = ("bakta", "smorfinder", "gmscmapper")
# Names the merge already reads as another tool (bakta writes 'Pyrodigal'), or uses itself.
RESERVED = set(BUILTIN) | {"orfipy"} | {k for k, v in SOURCE_TO_TOOL.items() if k != v}
KEYS = {"name", "env", "command", "description", "batch", "cpus", "memory", "time",
        "timeout", "version_cmd"}
DEFAULTS = {"batch": 20, "cpus": 1, "memory": "4 GB", "time": "4h", "timeout": "1h"}
CATEGORIES = ("smORF prediction", "intergenic smORF", "coding potential", "spurious tag")


def pipeline_defaults(config):
    """key -> value for the simple `name = value` lines of nextflow.config's params block."""
    vals, inside = {}, False
    for line in Path(config).read_text().splitlines():
        if re.match(r"^\s*params\s*\{", line):
            inside = True
            continue
        if inside and re.match(r"^\}", line):
            break
        m = re.match(r"^\s*(\w+)\s*=\s*('[^']*'|\"[^\"]*\"|[^\s/]+)", line) if inside else None
        if m:
            vals[m.group(1)] = yaml.safe_load(m.group(2).replace('"', "'"))
    return vals


def run_params(run_dir):
    params = pipeline_defaults(PIPE / "nextflow.config")
    if (PIPE / "conf/site.config").exists():
        params.update(pipeline_defaults(PIPE / "conf/site.config"))
    params.update(yaml.safe_load((Path(run_dir) / "params.yaml").read_text()) or {})
    return params


def check_callers(callers):
    """-> list of problems, one plain sentence each."""
    if callers in (None, []):
        return []
    if not isinstance(callers, list):
        return ["extra_callers must be a list of tool entries"]
    errors, seen = [], set()
    for i, c in enumerate(callers, 1):
        where = f"extra_callers #{i}"
        if not isinstance(c, dict):
            errors.append(f"{where}: must be a map with name, env and command")
            continue
        name = str(c.get("name", ""))
        where = f"extra_callers '{name or i}'"
        if not re.fullmatch(r"[a-z][a-z0-9_]{0,30}", name):
            errors.append(f"{where}: name must be lowercase letters, digits or _, starting with a letter")
        elif name in RESERVED:
            errors.append(f"{where}: name is taken by a built-in tool ({', '.join(sorted(RESERVED))})")
        elif name in seen:
            errors.append(f"{where}: name used twice")
        seen.add(name)
        unknown = set(c) - KEYS
        if unknown:
            errors.append(f"{where}: unknown keys {', '.join(sorted(unknown))}")
        env = str(c.get("env", "none"))
        if not re.fullmatch(r"module:[\w./,+:=-]+|(conda|sif):/[\w./,+:=-]+|none", env):
            errors.append(f"{where}: env must be module:<name>, conda:/<env path>, sif:/<image> or none")
        elif env.startswith("conda:") and not Path(env[6:], "bin").is_dir():
            errors.append(f"{where}: conda env {env[6:]} not found (no bin/)")
        elif env.startswith("sif:") and not Path(env[4:]).is_file():
            errors.append(f"{where}: image {env[4:]} not found")
        cmd = str(c.get("command", ""))
        for ph in ("{fasta}", "{gff}"):
            if ph not in cmd:
                errors.append(f"{where}: command must use {ph}")
        if re.search(r"['\"]\{(fasta|gff|cpus|gid|tmp)\}['\"]", cmd):
            errors.append(f"{where}: do not quote placeholders; they are quoted for you")
        for key in ("batch", "cpus"):
            if key in c and not (isinstance(c[key], int) and c[key] >= 1):
                errors.append(f"{where}: {key} must be a whole number >= 1")
        if "memory" in c and not re.fullmatch(r"\d+(\.\d+)? ?(MB|GB)", str(c["memory"])):
            errors.append(f"{where}: memory like '4 GB'")
        for key in ("time", "timeout"):
            if key in c and not re.fullmatch(r"\d+[smhd]", str(c[key])):
                errors.append(f"{where}: {key} like '90m' or '2h'")
    return errors


def conda_version(binary, package):
    """Version of a conda package, from the conda-meta of the env holding binary."""
    try:
        hits = sorted(Path(binary).resolve().parents[1].glob(f"conda-meta/{package}-[0-9]*.json"))
    except IndexError:
        hits = []
    return hits[-1].name[len(package) + 1:].split("-")[0] if hits else "version unknown"


def antifam_release(hmm):
    """'release 8.0, 278 families' from the AntiFam download's version file."""
    try:
        text = (Path(hmm).parent / "version").read_text()
    except OSError:
        return str(hmm)
    rel = re.search(r"release:\s*(\S+)", text)
    fams = re.search(r"families:\s*(\d+)", text)
    return f"AntiFam release {rel.group(1) if rel else '?'}, {fams.group(1) if fams else '?'} families"


def git_version(directory):
    try:
        out = subprocess.run(["git", "-C", str(directory), "log", "-1", "--format=%h %cs"],
                             capture_output=True, text=True, timeout=20).stdout.strip()
    except (OSError, subprocess.TimeoutExpired):
        out = ""
    return f"git {out}" if out else "local copy"


def tool_rows(p):
    rows = [
        {"category": "smORF prediction", "stage": "call", "tool": "bakta", "type": "built-in", "software": p.get("bakta_module"),
         "database": p.get("bakta_db"), "runs": f"per genome, {p.get('bakta_batch')} per task",
         "resources": "4 cores, 16 GB, ~20 min/genome",
         "finds": "all genes (Pyrodigal) plus Bakta's own small-ORF search"},
        {"category": "smORF prediction", "stage": "call", "tool": "smorfinder", "type": "built-in", "software": p.get("smorfinder_module"),
         "database": "built-in models", "runs": f"per genome, {p.get('smorf_batch')} per task",
         "resources": "1 core, 4 GB", "finds": "small ORFs scored by a deep-learning model"},
        {"category": "smORF prediction", "stage": "call", "tool": "gmscmapper", "type": "built-in", "software": p.get("gmsc_module"),
         "database": p.get("gmsc_db"), "runs": "per chunk of 500 genomes",
         "resources": f"16 cores, {p.get('gmsc_memory')}, {p.get('gmsc_queue')} queue",
         "finds": "ORFs matching the Global Microbial smORF Catalogue (high-quality hits)"},
    ]
    for c in p.get("extra_callers") or []:
        c = {**DEFAULTS, **c}
        rows.append({
            "category": "smORF prediction", "stage": "call", "tool": c["name"], "type": "extra", "software": c.get("env", "none"),
            "database": "-", "runs": f"per genome, {c['batch']} per task",
            "resources": f"{c['cpus']} core{'s' if c['cpus'] != 1 else ''}, {c['memory']}, {c['timeout']}/genome",
            "finds": c.get("description", "") or c["command"]})
    orfipy = p.get("orfipy_bin", "")
    rows += [
        {"category": "intergenic smORF", "stage": "combine", "tool": "orfipy", "type": "built-in",
         "software": f"orfipy {conda_version(orfipy, 'orfipy')}", "database": "-",
         "runs": "per chunk", "resources": "2 cores, 4 GB",
         "finds": f"intergenic smORFs: every {p.get('orf_min')}-{p.get('orf_max')} nt ORF "
                  f"({p.get('orf_starts')} starts) in DNA with no annotated gene (orfipy_novel)"},
        {"category": "spurious tag", "stage": "flags", "tool": "antifam", "type": "built-in",
         "software": p.get("hmmer_module"),
         "database": antifam_release(p["antifam_hmm"]) if p.get("antifam_hmm") else "not set (conf/site.config)",
         "runs": "per chunk, each distinct protein once", "resources": "4 cores, 4 GB",
         "finds": "flags proteins matching known spurious-ORF families (flag, not a filter)"},
        {"category": "coding potential", "stage": "flags", "tool": "popcorn", "type": "built-in",
         "software": f"Popcorn ({git_version(p['popcorn_dir']) if p.get('popcorn_dir') else 'not set (conf/site.config)'})",
         "database": "built-in model", "runs": "per chunk, each distinct DNA sequence once",
         "resources": "4 cores, 8 GB",
         "finds": "coding probability of each smORF's DNA (flag, not a filter)"},
    ]
    return rows


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("run_dir")
    ap.add_argument("--tsv")
    ap.add_argument("--shell", metavar="NAME")
    args = ap.parse_args()

    p = run_params(args.run_dir)
    errors = check_callers(p.get("extra_callers"))
    if errors:
        print("extra_callers problems:\n" + "\n".join(f"- {e}" for e in errors), file=sys.stderr)
        sys.exit(2)

    if args.shell:
        hits = [c for c in p.get("extra_callers") or [] if c["name"] == args.shell]
        if not hits:
            sys.exit(f"no extra caller named {args.shell} in params.yaml")
        c = {**DEFAULTS, **hits[0]}
        for key in ("name", "env", "command", "cpus", "timeout", "version_cmd"):
            print(f"{key}={shlex.quote(str(c.get(key, '') or ''))}")
        print(f"python_bin={shlex.quote(str(p.get('python_bin', 'python3')))}")
        return

    rows = sorted(tool_rows(p), key=lambda r: CATEGORIES.index(r["category"]))
    cols = ["category", "stage", "tool", "type", "software", "database", "runs", "resources", "finds"]
    print("| " + " | ".join(cols) + " |")
    print("|" + "---|" * len(cols))
    for r in rows:
        print("| " + " | ".join(str(r[c]).replace("|", "\\|") for c in cols) + " |")
    if args.tsv:
        with open(args.tsv, "w") as fh:
            fh.write("\t".join(cols) + "\n")
            for r in rows:
                fh.write("\t".join(str(r[c]) for c in cols) + "\n")


if __name__ == "__main__":
    main()

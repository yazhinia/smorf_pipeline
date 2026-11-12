#!/usr/bin/env python3
"""Record a checkpoint decision, so developers can see how users steer the pipeline.

    log_decision.py --run-dir RUN --stage merge --question "Continue to combine?" \
        --choice continue [--options "continue|adjust|stop"] [--comment "..."] \
        [--params-diff "chunk_size: 500 -> 100"] [--feedback-dir DIR]

Appends one row to <run>/run_log/decisions.tsv, the run's only decision log. The columns a
person reads come first (time, stage, question, choice, comment, params change), the ones
developers use to compare runs last; it opens as a table in a spreadsheet, or in a terminal with
`column -t -s $'\t' decisions.tsv | less -S`. With --feedback-dir (or params.yaml feedback_dir),
also appends the row to <feedback-dir>/decisions_all.tsv, the shared file developers read
across runs. An existing file keeps its own column order.
Uses only the standard library, so it runs with any python3.
"""
import argparse
import getpass
import json
import os
import re
from datetime import datetime
from pathlib import Path

COLUMNS = ["timestamp", "stage", "question", "choice", "comment", "params_diff",
           "options", "user", "run_dir", "params_md5", "pipeline_git"]


def clean(s):
    return re.sub(r"[\t\r\n]+", " ", s or "").strip()


def append_row(path, row):
    """Append in the file's own column order: logs written before the reorder stay aligned."""
    cols = COLUMNS
    if path.exists() and path.stat().st_size:
        with open(path) as fh:
            cols = fh.readline().rstrip("\n").split("\t")
    with open(path, "a") as fh:
        if cols is COLUMNS:
            fh.write("\t".join(COLUMNS) + "\n")
        fh.write("\t".join(clean(row.get(c, "")) for c in cols) + "\n")


def feedback_dir_from_params(run_dir):
    p = run_dir / "params.yaml"
    if p.exists():
        m = re.search(r"^feedback_dir:\s*['\"]?([^'\"\n#]+)", p.read_text(), re.M)
        if m and m.group(1).strip() not in ("null", "~", ""):
            return m.group(1).strip()
    return None


def main():
    ap = argparse.ArgumentParser(description=__doc__,
                                 formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--run-dir", required=True, type=Path)
    ap.add_argument("--stage", required=True)
    ap.add_argument("--question", required=True)
    ap.add_argument("--choice", required=True)
    ap.add_argument("--options", default="")
    ap.add_argument("--comment", default="")
    ap.add_argument("--params-diff", default="")
    ap.add_argument("--feedback-dir")
    args = ap.parse_args()

    run_dir = args.run_dir.resolve()
    state = {}
    if (run_dir / "run_state.json").exists():
        state = json.loads((run_dir / "run_state.json").read_text())
    row = {
        "timestamp": datetime.now().strftime("%Y-%m-%d %H:%M:%S"),
        "user": getpass.getuser(),
        "run_dir": str(run_dir),
        "stage": args.stage,
        "question": args.question,
        "options": args.options,
        "choice": args.choice,
        "comment": args.comment,
        "params_diff": args.params_diff,
        "params_md5": state.get("params_md5", ""),
        "pipeline_git": state.get("pipeline_git", ""),
    }
    log_dir = run_dir / "run_log"
    log_dir.mkdir(parents=True, exist_ok=True)
    append_row(log_dir / "decisions.tsv", row)

    shared = args.feedback_dir or feedback_dir_from_params(run_dir)
    if shared:
        try:
            Path(shared).mkdir(parents=True, exist_ok=True)
            append_row(Path(shared) / "decisions_all.tsv", row)
        except OSError as exc:
            print(f"WARN: could not write shared feedback to {shared}: {exc}")
    print(f"logged: {args.stage} -> {args.choice}")


if __name__ == "__main__":
    main()

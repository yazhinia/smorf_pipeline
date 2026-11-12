#!/bin/bash
# Try one extra caller from a run's params.yaml on one genome, before launching `call`.
# Runs the same script the pipeline task runs (run_caller.sh). Must run inside an LSF job.
#
#   test_caller.sh <run_dir> <name> [genome.fna]
#
# Default genome: the run's first prepared genome (needs `prepare` done), whose contigs are
# named <gid>_cNNN as the pipeline expects. Output: <run_dir>/caller_tests/<name>/.
set -euo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
die() { echo "test_caller.sh: $*" >&2; exit 2; }
[[ $# -ge 2 ]] || { sed -n '2,8p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
[[ -n "${LSB_JOBID:-}" ]] || die "run inside an LSF job: bsub -I -q normal -n 1 -M 4000 -R 'select[mem>4000] rusage[mem=4000]' $0 $*"

RUN_DIR="$(cd "$1" && pwd -P)"
NAME="$2"
cfgs=("${BIN}/../nextflow.config")
[[ -f "${BIN}/../conf/site.config" ]] && cfgs=("${BIN}/../conf/site.config" "${cfgs[@]}")
py="$(sed -n "/^ *python_bin *=/{s/^ *python_bin *= *'\([^']*\)'.*/\1/p;q}" "${cfgs[@]}")"
vars="$("${py}" "${BIN}/list_tools.py" "${RUN_DIR}" --shell "${NAME}")" || die "fix params.yaml first"
eval "${vars}"

fa="${3:-$(find "${RUN_DIR}/work" -path '*/chunk_*.genomes/*.fna' -print -quit 2>/dev/null || true)}"
[[ -n "${fa}" && -f "${fa}" ]] || die "no prepared genome; run the prepare stage first or give a FASTA"
fa="$(readlink -f "${fa}")"

out="${RUN_DIR}/caller_tests/${NAME}"
rm -rf "${out}"; mkdir -p "${out}"
cd "${out}"
printf '%s\n' "${command}" > command.txt
"${BIN}/run_caller.sh" --name "${name}" --env "${env}" --command-file command.txt \
  --out-dir result --cpus "${cpus}" --timeout "${timeout}" --python "${python_bin}" \
  ${version_cmd:+--version-cmd "${version_cmd}"} --keep "${fa}"

gid="$(basename "${fa}" .fna)"
echo "== genome: ${fa}"
echo "== status";  cat result/status.tsv
[[ -f result/version.txt ]] && { echo "== version"; cat result/version.txt; }
if [[ -f "result/${gid}.gff" ]]; then
  awk -F'\t' '!/^#/ {n++; if (($5 - $4 + 1) / 3 - 1 <= 100) s++}
    END {printf "== CDS: %d, of which <= 100 aa: %d\n", n, s}' "result/${gid}.gff"
  echo "== first rows"; grep -v '^#' "result/${gid}.gff" | head -3
fi
echo "== tool log (last lines)"; tail -5 "w_${name}_${gid}/tool.log" 2>/dev/null || true
echo "== files in ${out}"

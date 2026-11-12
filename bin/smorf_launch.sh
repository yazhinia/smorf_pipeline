#!/bin/bash
# Submit the Nextflow head job for one stage of a run, with bsub (never on the head node).
#
#   smorf_launch.sh <run_dir> <stop_after> [extra nextflow args...]
#
# <run_dir> is the run's own folder on /lustre or /data/pam, one per genome set; several runs
# can go at once in different folders. It must hold params.yaml (the run's options);
# Nextflow's cache (.nextflow/), work/ and results/ all live in it, so every relaunch with
# -resume reuses the tasks that already finished.
# Runs of <= small_run genomes (50) start sooner: the head job goes to the normal queue and runs
# the light Python steps itself (profile 'small'). SMORF_HEAD_QUEUE overrides the queue.
# Writes <run_dir>/run_state.json with the job id, stage and options hash.
set -euo pipefail

usage() { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
[[ $# -ge 2 ]] || usage
# -P: physical paths, so a /data/pam/... link into /lustre is seen as /lustre
RUN_DIR="$(cd "$1" && pwd -P)"
STOP_AFTER="$2"
shift 2

PIPELINE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd -P)"
NXF_MODULE="${SMORF_NXF_MODULE:-cellgen/nextflow/25.10.4}"
HEAD_MEM="${SMORF_HEAD_MEM:-8000}"

case "${STOP_AFTER}" in
  prepare|call|merge|combine|flags) ;;
  *) echo "stop_after must be prepare|call|merge|combine|flags" >&2; exit 2 ;;
esac
[[ -f "${RUN_DIR}/params.yaml" ]] || { echo "no params.yaml in ${RUN_DIR}" >&2; exit 2; }
[[ "${RUN_DIR}" == /lustre/* || "${RUN_DIR}" == /data/pam/* ]] || {
  echo "run dir ${RUN_DIR} must be on /lustre or /data/pam (PaM containers cannot see other paths)" >&2; exit 2; }

# Genomes in the run: from the prepare summary once prepare has run, else counted from the input.
genome_count() {
  local summary="${RUN_DIR}/results/qc/prepare_summary.tsv" input mode
  if [[ -s "${summary}" ]]; then
    awk -F'\t' '$1 == "genomes_ok" { print $2 }' "${summary}"; return
  fi
  input="$(sed -n 's/^input:[[:space:]]*//p' "${RUN_DIR}/params.yaml" | sed 's/[[:space:]]*#.*//; s/["'\'']//g')"
  mode="$(sed -n 's/^split_mode:[[:space:]]*\([a-z]*\).*/\1/p' "${RUN_DIR}/params.yaml")"
  [[ "${input}" == /* ]] || input="${RUN_DIR}/${input}"
  if [[ -d "${input}" ]]; then
    find "${input}" -type f -regextype posix-extended -iregex '.*\.(fa|fna|fasta)(\.gz)?' | wc -l
  elif [[ "${input}" =~ \.(tsv|csv)$ ]]; then
    echo $(( $(grep -c . "${input}") - 1 ))
  elif [[ "${mode:-file}" == file ]]; then
    echo 1
  else
    echo 999999      # one FASTA split into genomes by record: counted by prepare
  fi
}
small_run="$(sed -n 's/^small_run:[[:space:]]*\([0-9]*\).*/\1/p' "${RUN_DIR}/params.yaml")"
[[ -n "${small_run}" ]] || small_run="$(sed -n 's/^ *small_run *= *\([0-9]*\).*/\1/p' "${PIPELINE_DIR}/nextflow.config")"
n_genomes="$(genome_count)"
if (( ${n_genomes:-999999} <= ${small_run:-50} )); then
  profile="sanger,small"; HEAD_QUEUE="${SMORF_HEAD_QUEUE:-normal}"; xmx=3g
else
  profile="sanger"; HEAD_QUEUE="${SMORF_HEAD_QUEUE:-long}"; xmx=6g
fi

mkdir -p "${RUN_DIR}/run_log"
stamp="$(date '+%Y%m%dT%H%M%S')"
log="${RUN_DIR}/run_log/head_${STOP_AFTER}_${stamp}.log"
resume=""
[[ -d "${RUN_DIR}/.nextflow" ]] && resume="-resume"

cmd="set +u; source /etc/profile.d/modules.sh; module use /software/modules; \
module load ${NXF_MODULE}; set -u; cd '${RUN_DIR}'; \
export NXF_ANSI_LOG=false NXF_OPTS='-Xms1g -Xmx${xmx}'; \
nextflow run '${PIPELINE_DIR}/main.nf' -profile ${profile} -params-file params.yaml \
  -w '${RUN_DIR}/work' --outdir '${RUN_DIR}/results' --stop_after ${STOP_AFTER} ${resume} $*"

out="$(bsub -q "${HEAD_QUEUE}" -n 2 -M "${HEAD_MEM}" \
  -R "select[mem>${HEAD_MEM}] rusage[mem=${HEAD_MEM}]" \
  -J "smorf_${STOP_AFTER}_$(basename "${RUN_DIR}")" -o "${log}" -e "${log}" \
  bash -c "${cmd}")"
job_id="$(sed -n 's/^Job <\([0-9]*\)>.*/\1/p' <<<"${out}")"
[[ -n "${job_id}" ]] || { echo "bsub failed: ${out}" >&2; exit 1; }

git_sha="$(git -C "${PIPELINE_DIR}" rev-parse --short HEAD 2>/dev/null || echo none)"
params_md5="$(md5sum "${RUN_DIR}/params.yaml" | cut -c1-12)"
cat > "${RUN_DIR}/run_state.json" <<EOF
{
  "run_dir": "${RUN_DIR}",
  "pipeline_dir": "${PIPELINE_DIR}",
  "pipeline_git": "${git_sha}",
  "stage": "${STOP_AFTER}",
  "genomes": "${n_genomes}",
  "profile": "${profile}",
  "head_queue": "${HEAD_QUEUE}",
  "job_id": "${job_id}",
  "head_log": "${log}",
  "params_md5": "${params_md5}",
  "resume": "${resume:-no}",
  "submitted": "${stamp}"
}
EOF
echo "submitted job ${job_id} (${HEAD_QUEUE} queue, profile ${profile}): stop_after=${STOP_AFTER} ${resume} log=${log}"

#!/bin/bash
# Run one extra gene caller (params.yaml: extra_callers) on genome FASTAs. Used by the
# CUSTOM_CALLER task and by test_caller.sh, so a test runs exactly what the pipeline runs.
#
#   run_caller.sh --name NAME --env ENV --command-file FILE --out-dir DIR \
#       [--cpus 1] [--timeout 1h] [--python python3] [--version-cmd CMD] [--keep] FASTA...
#
# ENV: module:<m1>[,<m2>] (PaM modules) | conda:<env prefix> | sif:<image.sif> | none
# FILE holds the tool's command; {fasta} {gff} {cpus} {gid} {tmp} are filled in per genome.
# Writes DIR/<gid>.gff (cleaned by custom_gff.py), DIR/status.tsv and, with --version-cmd,
# DIR/version.txt. A genome that fails is a status row; the run goes on.
set -uo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
die() { echo "run_caller.sh: $*" >&2; exit 2; }

name="" env="none" cmd_file="" out_dir="" cpus=1 limit=1h python=python3 version_cmd="" keep=""
while [[ $# -gt 0 && "$1" == --* ]]; do
  case "$1" in
    --name) name="$2"; shift 2 ;;
    --env) env="$2"; shift 2 ;;
    --command-file) cmd_file="$2"; shift 2 ;;
    --out-dir) out_dir="$2"; shift 2 ;;
    --cpus) cpus="$2"; shift 2 ;;
    --timeout) limit="$2"; shift 2 ;;
    --python) python="$2"; shift 2 ;;
    --version-cmd) version_cmd="$2"; shift 2 ;;
    --keep) keep=1; shift ;;
    *) die "unknown option $1" ;;
  esac
done
[[ -n "${name}" && -f "${cmd_file}" && -n "${out_dir}" && $# -gt 0 ]] \
  || die "needs --name, --command-file, --out-dir and at least one FASTA"
command -v "${python}" >/dev/null || die "python '${python}' not found (params python_bin)"

# Tool environment. Module or image problems stop the task: they would fail every genome.
wrap=()
case "${env}" in
  module:*)
    IFS=',' read -ra mods <<<"${env#module:}"
    source "${BIN}/pam_modules.sh" "${mods[@]}" ;;
  conda:*)
    prefix="${env#conda:}"
    [[ -d "${prefix}/bin" ]] || { echo "conda env ${prefix} has no bin/" >&2; exit 3; }
    export PATH="${prefix}/bin:${PATH}" CONDA_PREFIX="${prefix}" ;;
  sif:*)
    img="${env#sif:}"
    [[ -f "${img}" ]] || { echo "image ${img} not found" >&2; exit 3; }
    source "${BIN}/pam_modules.sh" "${SMORF_SINGULARITY_MODULE:-singularityce-4.1.0}"
    wrap=(singularity exec -B /lustre,/data/pam "${img}") ;;
  none|"") ;;
  *) die "env must be module:..., conda:..., sif:... or none; got '${env}'" ;;
esac

mkdir -p "${out_dir}"
out_dir="$(cd "${out_dir}" && pwd -P)"
printf 'genome_id\ttool\tstatus\treason\tn_records\n' > "${out_dir}/status.tsv"

cmd="${PWD}/caller_cmd.sh"
sed -e 's/{fasta}/"$FASTA"/g' -e 's/{gff}/"$GFF"/g' -e 's/{cpus}/"$CPUS"/g' \
    -e 's/{gid}/"$GID"/g' -e 's/{tmp}/"$TMPDIR"/g' "${cmd_file}" > "${cmd}"

if [[ -n "${version_cmd}" ]]; then
  v="$(${wrap[@]+"${wrap[@]}"} bash -c "${version_cmd}" 2>&1 | grep -m1 . || true)"
  printf '%s\t%s\n' "${name}" "${v:-unknown}" > "${out_dir}/version.txt"
fi

for fa in "$@"; do
  fa="$(readlink -f "${fa}")"
  gid="$(basename "${fa}" .fna)"
  w="${PWD}/w_${name}_${gid}"
  rm -rf "${w}"; mkdir -p "${w}/tmp"
  rc=0
  ( cd "${w}" && export FASTA="${fa}" GFF="${w}/out.gff" GID="${gid}" CPUS="${cpus}" \
        TMPDIR="${w}/tmp" SINGULARITYENV_TMPDIR="${w}/tmp" \
      && timeout "${limit}" ${wrap[@]+"${wrap[@]}"} bash "${cmd}" ) > "${w}/tool.log" 2>&1 || rc=$?
  "${python}" "${BIN}/custom_gff.py" --tool "${name}" --genome "${fa}" --raw "${w}/out.gff" \
      --rc "${rc}" --out "${out_dir}/${gid}.gff" --status "${out_dir}/status.tsv"
  [[ -f "${out_dir}/${gid}.gff" ]] || tail -5 "${w}/tool.log" >&2 || true
  [[ -n "${keep}" ]] || rm -rf "${w}"
done

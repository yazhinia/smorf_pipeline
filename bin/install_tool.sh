#!/bin/bash
# Install an extra gene caller into a shared tools folder, never into the pipeline copy.
# Must run inside an LSF job (bsub -I ...): conda solves and image pulls are too heavy for
# the head node.
#
#   install_tool.sh <tools_dir> <name> conda <package spec>...   e.g. pyrodigal=3.6.3 (bioconda)
#   install_tool.sh <tools_dir> <name> pip <package spec>...     Python env + pip packages
#   install_tool.sh <tools_dir> <name> sif <image uri>           e.g. docker://quay.io/biocontainers/...
#
# Installs into <tools_dir>/<name>/ and prints the `env:` value for params.yaml extra_callers.
# <tools_dir> must be on /lustre or /data/pam so every farm node (and container) can read it.
set -euo pipefail

BIN="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd -P)"
die() { echo "install_tool.sh: $*" >&2; exit 2; }
[[ $# -ge 4 ]] || { sed -n '2,12p' "$0" | sed 's/^# \{0,1\}//'; exit 2; }
[[ -n "${LSB_JOBID:-}" ]] || die "run inside an LSF job: bsub -I -q normal -n 2 -M 8000 -R 'select[mem>8000] rusage[mem=8000]' $0 $*"

mkdir -p "$1"
TOOLS_DIR="$(cd "$1" && pwd -P)"
NAME="$2"
METHOD="$3"
shift 3
[[ "${TOOLS_DIR}" == /lustre/* || "${TOOLS_DIR}" == /data/pam/* ]] || die "tools dir must be on /lustre or /data/pam"
[[ "${NAME}" =~ ^[a-z][a-z0-9_]{0,30}$ ]] || die "name must be lowercase letters, digits or _"
DEST="${TOOLS_DIR}/${NAME}"
[[ ! -e "${DEST}" ]] || die "${DEST} exists; remove it or pick another name"
mkdir -p "${DEST}"
export TMPDIR="${DEST}/.tmp"
mkdir -p "${TMPDIR}"

conda_bin="${SMORF_CONDA:-$(command -v mamba || command -v conda || echo /software/conda/miniforge/bin/mamba)}"
export CONDA_PKGS_DIRS="${TOOLS_DIR}/.conda_pkgs"   # off the /nfs home quota

case "${METHOD}" in
  conda)
    "${conda_bin}" create -y -p "${DEST}/env" --override-channels -c conda-forge -c bioconda "$@"
    "${conda_bin}" list -p "${DEST}/env" > "${DEST}/packages.txt"
    env_line="conda:${DEST}/env" ;;
  pip)
    "${conda_bin}" create -y -p "${DEST}/env" --override-channels -c conda-forge python=3.11 pip
    "${DEST}/env/bin/pip" install --no-cache-dir "$@"
    "${DEST}/env/bin/pip" freeze > "${DEST}/packages.txt"
    env_line="conda:${DEST}/env" ;;
  sif)
    [[ $# -eq 1 ]] || die "sif takes one image uri"
    source "${BIN}/pam_modules.sh" "${SMORF_SINGULARITY_MODULE:-singularityce-4.1.0}"
    export SINGULARITY_CACHEDIR="${TOOLS_DIR}/.sif_cache"
    singularity pull "${DEST}/${NAME}.sif" "$1"
    echo "$1" > "${DEST}/packages.txt"
    env_line="sif:${DEST}/${NAME}.sif" ;;
  *) rm -rf "${DEST}"; die "method must be conda, pip or sif" ;;
esac
rm -rf "${TMPDIR}"

printf 'date\tuser\tmethod\tspec\tenv\n%s\t%s\t%s\t%s\t%s\n' \
  "$(date '+%Y-%m-%dT%H:%M:%S')" "${USER:-unknown}" "${METHOD}" "$*" "${env_line}" > "${DEST}/INSTALL.tsv"
echo "installed ${NAME}: env: ${env_line}"

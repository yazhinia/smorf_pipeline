# Sourced, not executed:  source pam_modules.sh <module> [<module> ...]
#
# LSF runs task scripts in a plain non-login shell, so neither the `module` function nor the
# PaM module path exists yet. Same setup as bakta_worker.sh / smorfinder_worker.sh.
# A module that fails to load must stop the task here, not surface later as an empty output.
set +u
if ! command -v module >/dev/null 2>&1; then
  source /etc/profile.d/modules.sh
fi
module use "${SMORF_MODULE_DIR:-/software/modules}"
module load "${SMORF_ENV_MODULE:-PaM/environment}"
for _m in "$@"; do
  module load "${_m}" || { echo "module load ${_m} failed" >&2; exit 3; }
done
unset _m
set -u

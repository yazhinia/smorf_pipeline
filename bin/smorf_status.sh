#!/bin/bash
# One-screen progress report for a run started by smorf_launch.sh.
#
#   smorf_status.sh <run_dir>
#
# Prints the head job state, task counts per process from the Nextflow trace, failed tasks,
# and which stage summaries exist. Exit code: 0 finished OK, 1 failed, 3 still running.
set -uo pipefail

RUN_DIR="$(cd "${1:?usage: smorf_status.sh <run_dir>}" && pwd -P)"
state="${RUN_DIR}/run_state.json"
[[ -f "${state}" ]] || { echo "no run_state.json in ${RUN_DIR}: nothing launched yet"; exit 1; }

field() { sed -n "s/.*\"$1\": \"\([^\"]*\)\".*/\1/p" "${state}"; }
job_id="$(field job_id)"; stage="$(field stage)"; log="$(field head_log)"

job_state="$(bjobs -noheader -o stat "${job_id}" 2>/dev/null | head -1)"
[[ -n "${job_state}" ]] || job_state="$(bhist -n 0 -l "${job_id}" 2>/dev/null | grep -o 'Done successfully\|Exited' | tail -1)"
[[ -n "${job_state}" ]] || job_state=UNKNOWN
echo "run: ${RUN_DIR}"
echo "stage target: ${stage}   head job: ${job_id} (${job_state})"

trace="${RUN_DIR}/results/run_log/trace.tsv"
if [[ -s "${trace}" ]]; then
  echo
  echo "tasks per process (status counts):"
  awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    {p=$c["process"]; sub(/.*:/,"",p); n[p" "$c["status"]]++; procs[p]=1}
    END{for(p in procs){line=sprintf("  %-16s",p);
        for(k in n){split(k,a," "); if(a[1]==p) line=line sprintf(" %s=%d",a[2],n[k])}
        print line}}' "${trace}" | sort
  failed="$(awk -F'\t' 'NR==1{for(i=1;i<=NF;i++) c[$i]=i; next}
    $c["status"]=="FAILED"{print "  " $c["process"] " " $c["tag"] " exit=" $c["exit"] " " $c["workdir"]}' "${trace}")"
  if [[ -n "${failed}" ]]; then echo; echo "failed tasks (see .command.err in the work dir):"; echo "${failed}" | head -10; fi
fi

echo
echo "stage summaries:"
for s in prepare call merge combine flags; do
  md="${RUN_DIR}/results/qc/${s}_summary.md"
  if [[ -f "${md}" ]]; then echo "  ${s}: $(head -1 "${md}" | sed 's/^## //')"; else echo "  ${s}: -"; fi
done

if [[ -f "${log}" ]] && grep -q 'ERROR ~' "${log}"; then
  echo; echo "head log errors:"; grep -A5 'ERROR ~' "${log}" | head -12
fi

case "${job_state}" in
  DONE|"Done successfully") exit 0 ;;
  EXIT|Exited) exit 1 ;;
  *) exit 3 ;;
esac

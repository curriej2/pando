#!/usr/bin/env bash
# Submit the rung-2 silencing library (24_silencing_library.py) as a 30-task array, one tree cell each.
# Run ONLY after 23_validate_silencing.py has PASSED (CLAUDE.md, PROPOSAL (2026-10-01, rev. 1) §4b).
# Pure CPU, one core per task: partition cpu, never lesliec (CLAUDE.md "Do not list lesliec ...").
#   analyses/2026-09_simulator/src/25_submit_silencing_library.sh [--mem 4G] [--time 01:00:00] [--array 0-29] [--quick]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
MEM=4G; TIME=01:00:00; ARRAY=0-29; EXTRA=""
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --mem)   MEM="$2";   shift 2;;
    --time)  TIME="$2";  shift 2;;
    --array) ARRAY="$2"; shift 2;;
    --quick) EXTRA="--quick"; shift;;
    *) echo "unknown flag $1" >&2; exit 1;;
  esac
done
mkdir -p "$ROOT/logs"
sbatch -A lesliec -p cpu -c 1 --mem "$MEM" -t "$TIME" --array "$ARRAY" \
       -J silencing_library -o "$ROOT/logs/%x-%A_%a.out" --wrap \
       "cd $ROOT && /data1/choij10/justin/envs/pando/bin/python -u analyses/2026-09_simulator/src/24_silencing_library.py \$SLURM_ARRAY_TASK_ID $EXTRA"

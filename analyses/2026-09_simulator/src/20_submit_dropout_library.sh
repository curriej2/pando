#!/usr/bin/env bash
# Submit the rung-1 dropout library (19_dropout_library.py) as a 30-task array, one tree cell each.
# Run ONLY after 18_validate_dropout.py has PASSED (CLAUDE.md, PROPOSAL (2026-09-29), amendment 6).
# Pure CPU, one core per task: partition cpu, never lesliec (CLAUDE.md "Do not list lesliec ...").
#   analyses/2026-09_simulator/src/20_submit_dropout_library.sh [--mem 4G] [--time 01:00:00] [--array 0-29]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/../../.." && pwd)"
MEM=4G; TIME=01:00:00; ARRAY=0-29
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --mem)   MEM="$2";   shift 2;;
    --time)  TIME="$2";  shift 2;;
    --array) ARRAY="$2"; shift 2;;
    *) echo "unknown flag $1" >&2; exit 1;;
  esac
done
mkdir -p "$ROOT/logs"
sbatch -A lesliec -p cpu -c 1 --mem "$MEM" -t "$TIME" --array "$ARRAY" \
       -J dropout_library -o "$ROOT/logs/%x-%A_%a.out" --wrap \
       "cd $ROOT && /data1/choij10/justin/envs/pando/bin/python -u analyses/2026-09_simulator/src/19_dropout_library.py \$SLURM_ARRAY_TASK_ID"

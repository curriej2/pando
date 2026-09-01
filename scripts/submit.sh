#!/usr/bin/env bash
# Submit an analysis script as a Slurm batch job. Defaults are deliberately MODEST
# (16G, 4h) so jobs backfill instead of queueing: raise them from a measured MaxRSS
# (sacct -j <id> -o MaxRSS), not from a guess. See CLAUDE.md "Right-size the request".
# run in the VS Code tunnel -- a crash there kills the session AND the computation.
#
#   scripts/submit.sh analyses/2026-08_park-compatibility/src/14_compat_sparse.py Mouse1 7200
#   scripts/submit.sh --mem 512G --time 2-00:00:00 <script> [args...]
set -euo pipefail
ROOT="$(cd "$(dirname "${BASH_SOURCE[0]}")/.." && pwd)"
MEM=16G; TIME=4:00:00; CPUS=4; DEP=""   # modest by default -- raise only from a measured MaxRSS
while [[ "${1:-}" == --* ]]; do
  case "$1" in
    --mem)  MEM="$2";  shift 2;;
    --time) TIME="$2"; shift 2;;
    --cpus) CPUS="$2"; shift 2;;
    --dep)  DEP="--dependency=afterok:$2"; shift 2;;
    *) echo "unknown flag $1" >&2; exit 1;;
  esac
done
SCRIPT="${1:?usage: submit.sh [--mem M] [--time T] [--cpus N] <script.py> [args...]}"; shift
NAME="$(basename "$SCRIPT" .py)${1:+_$1}"
mkdir -p "$ROOT/logs"
sbatch $DEP -A lesliec -p lesliec,cpu -c "$CPUS" --mem "$MEM" -t "$TIME" \
       -J "$NAME" -o "$ROOT/logs/%x-%j.out" --wrap \
       "cd $ROOT && /data1/choij10/justin/envs/pando/bin/python -u $SCRIPT $*"

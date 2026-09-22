#!/usr/bin/env bash
# Submit 04_tree_library.py as a Slurm ARRAY.  scripts/submit.sh has no --array
# support, so this follows the precedent of 47_submit_B1000.sh in the
# park-compatibility analysis: a dedicated submitter for array work.
#
# Sizing, from 04's --dry-run: 554 chunks / 7,400 trees / 352 core-h, balanced
# across 120 tasks at min 1.82 h, median 2.57 h, max 4.91 h.  12 h walltime is
# ~2.4x the heaviest task.  Peak RSS measured <= 0.18 GB anywhere in this
# analysis, so 8 G is already generous.
#
# ⚠ -p cpu ALONE.  A wide array of small CPU tasks must never list lesliec:
# those four nodes are the lab's only GPU nodes (CLAUDE.md, 2026-09-03).
set -euo pipefail
ROOT=/data1/choij10/justin/pando
NTASKS=${1:-120}
REPS=${2:-20}
mkdir -p "$ROOT/logs"
sbatch -A lesliec -p cpu -c 1 --mem 8G -t 12:00:00 \
       --array=0-$((NTASKS-1)) -J tree_library \
       -o "$ROOT/logs/%x-%A_%a.out" --wrap \
       "cd $ROOT && /data1/choij10/justin/envs/pando/bin/python -u \
        analyses/2026-09_simulator/src/04_tree_library.py \
        --task \$SLURM_ARRAY_TASK_ID --ntasks $NTASKS --reps $REPS"

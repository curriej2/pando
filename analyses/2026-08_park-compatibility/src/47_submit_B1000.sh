#!/usr/bin/env bash
# =============================================================================
#  B = 1,000 permutations for the whole event-catalogue programme
# =============================================================================
#  Submits, for every arm:
#    * 42_clonewide.py  --nperm 1000                 one job   (~2 min, no split)
#    * 40_event_catalogue.py --permpart i/20         20 array tasks
#    * 43_soft_events.py --maxd 4 --permpart i/20    20 array tasks
#    * 43_soft_events.py --maxd 6 --permpart i/20    20 array tasks
#  then ONE dependent job that pools the parts (46_perm_merge.py, which refuses
#  an incomplete or duplicated set) and re-runs each observed scan once with
#  --nullfile, attaching q-values and the global p.
#
#  Sizing comes from a COMPLETED job's sacct, not a guess: peak RSS across every
#  40/43 run to date is 935 MB (43 Initial), so 8 G is ~8x headroom.  Walltimes
#  are 4-6x the measured per-scan cost x 50 permutations per task.  One core per
#  task: the inner loop is bincount/exp, not BLAS, so extra cores buy nothing and
#  a 1-core 8 G job backfills into gaps a fat one cannot.
#
#  ⚠ PARTITION: `cpu` ONLY -- deliberately not the usual `-p lesliec,cpu` pair.
#  300 one-core tasks listed against both landed 194 CPUs on lesliec, 76% of the
#  lab's four private nodes, held by one user.  Those four nodes are also the
#  lab's ONLY GPU nodes, so pure-CPU work parked on them can block a labmate's
#  A100 job on CPUs while the GPUs sit idle -- the worst way to use a private
#  allocation.  The general partition has 239 nodes and ~9,700 idle CPUs, 30x
#  what this needs, and a 1-core 8 G task is exactly the shape that backfills
#  there.  Keep lesliec for work that needs its GPUs or its ~1 TB nodes.
#
#  Usage: bash src/47_submit_B1000.sh [--dry] [--parts-only]
#         --parts-only skips the five 42_clonewide runs (2 min each, already done)
# =============================================================================
set -euo pipefail
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
B=1000; NP=20; SEED=20260903
ARMS="Mouse3 Mouse1 Mouse2 Initial Subclone"
DRY=0; SKIP42=0
for a in "$@"; do
  case "$a" in --dry) DRY=1;; --parts-only) SKIP42=1;;
     *) echo "unknown flag $a" >&2; exit 1;; esac
done
mkdir -p "$ROOT/logs" "$D/results/permparts"

sub() {  # sub <name> <time> <extra sbatch args...> -- <command>; echoes the job id
  local name="$1" tm="$2"; shift 2
  local extra=(); while [[ "$1" != "--" ]]; do extra+=("$1"); shift; done; shift
  if [[ $DRY == 1 ]]; then
    echo "[dry] sbatch --mem 8G -t $tm -J $name ${extra[*]-} --wrap '$*'" >&2
    echo 000000; return
  fi
  sbatch -A lesliec -p cpu -c 1 --mem 8G -t "$tm" -J "$name" \
       -o "$ROOT/logs/%x-%A_%a.out" "${extra[@]}" \
       --wrap "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $*" | awk '{print $NF}'
}

DEPS=()
for A in $ARMS; do
  # --- clone-wide layer: cheap enough to run B=1000 in one job
  [[ $SKIP42 == 0 ]] && sub "perm42_$A" 2:00:00 -- \
      "$D/src/42_clonewide.py $A --nperm $B --seed $SEED" >/dev/null

  # --- sub-clone hard catalogue
  DEPS+=( "$(sub "perm40_$A" 6:00:00 --array=1-$NP -- \
      "$D/src/40_event_catalogue.py $A --nperm $B --seed $SEED --permpart \$SLURM_ARRAY_TASK_ID/$NP")" )

  # --- sub-clone soft catalogue, at depth 4 and at full recorder depth 6
  for MD in 4 6; do
    T=6:00:00; [[ "$A" == Initial || "$A" == Subclone ]] && T=16:00:00
    DEPS+=( "$(sub "perm43d${MD}_$A" "$T" --array=1-$NP -- \
        "$D/src/43_soft_events.py $A --maxd $MD --nperm $B --seed $SEED --permpart \$SLURM_ARRAY_TASK_ID/$NP")" )
  done
done

# --- merge + re-run every observed scan once against the pooled null.
# afterany, not afterok: if one array task dies we want 46_perm_merge.py to say
# WHICH permutations are missing, not to have the collector silently cancelled.
DEP=$(IFS=:; echo "${DEPS[*]}")
[[ $DRY == 1 ]] && { echo "[dry] collector would depend on afterany:$DEP"; exit 0; }
sbatch -A lesliec -p cpu -c 1 --mem 16G -t 6:00:00 -J perm_collect \
     --dependency=afterany:"$DEP" -o "$ROOT/logs/%x-%j.out" \
     --wrap "cd $ROOT && export OMP_NUM_THREADS=1 && bash $D/src/48_collect_B1000.sh"
echo "submitted; collector depends on: $DEP"

#!/usr/bin/env bash
# =============================================================================
#  B = 1,000 at a LOW scan floor, with a clade-size-stratified null
# =============================================================================
#  Supersedes 47_submit_B1000.sh (kept as the record of the floor-10 run).
#  Two changes, both forced by the 2026-09-04 audit:
#
#  1. LOWER FLOOR.  At floor 10 the adaptive FDR<=5% step was degenerate in 13
#     of 15 configurations -- FDR there is 262-2,511x below target, so the
#     threshold simply snapped to an arbitrary constant. Hard goes to 2 nats,
#     soft to 4 (script 45 measured Lambda_soft's null mean at -1.74 sd 3.90, so
#     below ~4 nats the soft statistic carries no information worth calibrating).
#  2. STRATIFIED NULL.  A fixed Lambda floor is a clade-size filter in disguise:
#     a fully-missing m-cell clade at expected rate p caps near m*log(1/p) nats,
#     so a 9-cell Pre-TX clade cannot exceed ~12.4 however complete the loss.
#     Counts are now kept per size stratum and the threshold chosen per stratum.
#
#  Costs measured on the Mouse3 smoke test: the streaming histogram makes the
#  lower floor nearly free (1.0 s/scan at floor 2 vs 1.1 s at floor 10 for the
#  hard scan), because only the survivors of the same L >= LAM0 prefilter are
#  binned and nothing is concatenated.
#
#  Outputs are tagged _lam2 / _lam4 so the floor-10 catalogue is NOT overwritten.
#  Usage: bash src/49_submit_lowfloor.sh [--dry]
# =============================================================================
set -euo pipefail
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
B=1000; NP=20; SEED=20260903; LAM_HARD=2; LAM_SOFT=4
ARMS="Mouse3 Mouse1 Mouse2 Initial Subclone"
DRY=0; [[ "${1:-}" == "--dry" ]] && DRY=1
mkdir -p "$ROOT/logs" "$D/results/permparts"

sub() {  # sub <name> <mem> <time> <extra...> -- <command>;  echoes the job id
  local name="$1" mem="$2" tm="$3"; shift 3
  local extra=(); while [[ "$1" != "--" ]]; do extra+=("$1"); shift; done; shift
  if [[ $DRY == 1 ]]; then
    echo "[dry] sbatch --mem $mem -t $tm -J $name ${extra[*]-} --wrap '$*'" >&2
    echo 000000; return
  fi
  # -p cpu only: a wide array of one-core tasks must not sit on the lab's
  # four private GPU nodes (see root CLAUDE.md).
  sbatch -A lesliec -p cpu -c 1 --mem "$mem" -t "$tm" -J "$name" \
       -o "$ROOT/logs/%x-%A_%a.out" "${extra[@]}" \
       --wrap "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $*" | awk '{print $NF}'
}

DEPS=()
for A in $ARMS; do
  T=8:00:00; [[ "$A" == Initial || "$A" == Subclone ]] && T=24:00:00
  DEPS+=( "$(sub "lf40_$A" 8G "$T" --array=1-$NP -- \
      "$D/src/40_event_catalogue.py $A --lam $LAM_HARD --nperm $B --seed $SEED \
       --permpart \$SLURM_ARRAY_TASK_ID/$NP")" )
  for MD in 4 6; do
    DEPS+=( "$(sub "lf43d${MD}_$A" 12G "$T" --array=1-$NP -- \
        "$D/src/43_soft_events.py $A --maxd $MD --lam $LAM_SOFT --nperm $B --seed $SEED \
         --permpart \$SLURM_ARRAY_TASK_ID/$NP")" )
  done
done

DEP=$(IFS=:; echo "${DEPS[*]}")
[[ $DRY == 1 ]] && { echo "[dry] collector would depend on afterany:$DEP"; exit 0; }
# the collector runs the OBSERVED scans, which do hold per-candidate arrays at a
# low floor -- hence 64 G here against 8-12 G for the counts-only parts.
sbatch -A lesliec -p cpu -c 1 --mem 64G -t 12:00:00 -J lf_collect \
     --dependency=afterany:"$DEP" -o "$ROOT/logs/%x-%j.out" \
     --wrap "cd $ROOT && export OMP_NUM_THREADS=1 && bash $D/src/50_collect_lowfloor.sh"
echo "submitted; collector depends on: $DEP"

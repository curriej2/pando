#!/usr/bin/env bash
# Per-combo nulls, B=1000, 5 parts per arm.  Sized from measured runs (2026-09-07):
#   Pre-TX  103.3M slots, 36 s/scan, 4.87 GB peak, 1.09 GB per part on disk
#   Mouse3    5.6M slots, 0.6 s/scan, 0.34 GB peak
# 200 permutations per part => Pre-TX ~2 h, the mice minutes.  -p cpu only.
set -euo pipefail
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
B=1000; NP=5; SEED=20260903
DEPS=()
for A in Mouse3 Mouse1 Mouse2 Initial Subclone; do
  MEM=8G; T=4:00:00
  [[ "$A" == Subclone ]] && { MEM=16G; T=12:00:00; }
  [[ "$A" == Initial  ]] && { MEM=24G; T=12:00:00; }
  DEPS+=( "$(sbatch -A lesliec -p cpu -c 1 --mem $MEM -t $T -J "pc_$A" --array=1-$NP \
      -o "$ROOT/logs/%x-%A_%a.out" --wrap \
      "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $D/src/53_percombo.py $A \
       --nperm $B --seed $SEED --permpart \$SLURM_ARRAY_TASK_ID/$NP" | awk '{print $NF}')" )
done
DEP=$(IFS=:; echo "${DEPS[*]}")
sbatch -A lesliec -p cpu -c 1 --mem 48G -t 6:00:00 -J pc_merge \
   --dependency=afterany:"$DEP" -o "$ROOT/logs/%x-%j.out" --wrap \
   "cd $ROOT && export OMP_NUM_THREADS=1 && for A in Mouse3 Mouse1 Mouse2 Initial Subclone; do \
    $PY -u $D/src/54_percombo_merge.py \$A || echo FAILED \$A; done"
echo "submitted; merge depends on: $DEP"

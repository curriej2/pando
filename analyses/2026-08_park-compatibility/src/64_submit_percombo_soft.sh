#!/usr/bin/env bash
# Per-combo SOFT nulls at full recorder depth.  One arm per invocation.
#   accumulation : B permutations split over NP array tasks
#   calibration  : NCAL FRESH permutations (indices B..B+NCAL-1), one task
#   merge        : depends on both
# ⚠ -p cpu only: a wide array of small tasks must not park on the lab's 4 GPU nodes.
set -euo pipefail
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
A=${1:?arm}; B=${2:-1000}; NP=${3:-5}; NCAL=${4:-12}; MAXD=${5:-6}
# Sized from the measured Mouse3 run (77 B/slot in accumulate mode, plus codes,
# block index arrays and the per-cell copies) and the d6 slot counts:
#   Mouse3 6.85M · Mouse2 9.25M · Mouse1 19.06M · Subclone 61.73M · Pre-TX 105.88M
# ⚠ CALIB needs its own figure: it holds z for every draw (NCAL x NC float32 =
#   5.1 GB on Pre-TX at NCAL=12) but skips the observed pass.
case "$A" in
  Mouse3)   MEM=8G;  CMEM=8G;  T=3:00:00;  MMEM=16G ;;
  Mouse2)   MEM=8G;  CMEM=8G;  T=6:00:00;  MMEM=24G ;;
  Mouse1)   MEM=8G;  CMEM=12G; T=8:00:00;  MMEM=32G ;;
  Subclone) MEM=16G; CMEM=24G; T=24:00:00; MMEM=64G ;;
  Initial)  MEM=24G; CMEM=32G; T=24:00:00; MMEM=96G ;;
esac
ACC=$(sbatch --parsable -A lesliec -p cpu -c 1 --mem $MEM -t $T -J "ps_$A" --array=1-$NP \
  -o "$ROOT/logs/%x-%A_%a.out" --wrap \
  "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $D/src/62_percombo_soft.py $A \
   --nperm $B --maxd $MAXD --permpart \$SLURM_ARRAY_TASK_ID/$NP")
CAL=$(sbatch --parsable -A lesliec -p cpu -c 1 --mem $CMEM -t $T -J "pscal_$A" \
  -o "$ROOT/logs/%x-%j.out" --wrap \
  "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $D/src/62_percombo_soft.py $A \
   --nperm $B --maxd $MAXD --calib $NCAL")
MRG=$(sbatch --parsable -A lesliec -p cpu -c 1 --mem $MMEM -t 4:00:00 -J "psmrg_$A" \
  --dependency=afterok:$ACC:$CAL -o "$ROOT/logs/%x-%j.out" --wrap \
  "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $D/src/63_percombo_soft_merge.py $A \
   --maxd $MAXD --budget 2")
echo "$A: accum $ACC (array 1-$NP), calib $CAL, merge $MRG"

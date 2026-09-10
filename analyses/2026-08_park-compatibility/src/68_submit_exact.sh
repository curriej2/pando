#!/usr/bin/env bash
# Exact-moment analytic test on the STORED scan output -- no rescanning.
# Measured on Mouse3: 1.76 GB / 9 s at 6.85M slots = ~257 bytes per slot.
# d6 slot counts: Mouse3 6.85M · Mouse2 9.25M · Mouse1 19.06M · Subclone 61.73M · Pre-TX 105.88M
set -euo pipefail
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
EF=${1:-1}
for A in Mouse3 Mouse2 Mouse1 Subclone Initial; do
  case "$A" in
    Mouse3)   MEM=8G  ;;
    Mouse2)   MEM=12G ;;
    Mouse1)   MEM=20G ;;
    Subclone) MEM=48G ;;
    Initial)  MEM=72G ;;
  esac
  sbatch -A lesliec -p cpu -c 1 --mem $MEM -t 4:00:00 -J "ex_$A" \
    -o "$ROOT/logs/%x-%j.out" --wrap \
    "cd $ROOT && export OMP_NUM_THREADS=1 && $PY -u $D/src/67_exact_merge.py $A \
     --maxd 6 --efalse $EF --validate" | awk -v a=$A '{print a": "$NF}'
done

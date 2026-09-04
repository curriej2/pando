#!/usr/bin/env bash
# Pool the permutation parts and re-run each observed scan once against them.
# Run by the dependent job that 47_submit_B1000.sh submits; safe to rerun by hand.
# Deliberately NOT set -e: one arm failing must not stop the other fourteen, so
# failures are collected and reported at the end with a nonzero exit.
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
FAIL=()
cd "$ROOT"
for A in Mouse3 Mouse1 Mouse2 Initial Subclone; do
  echo "################################################## 40 $A"
  if $PY -u "$D/src/46_perm_merge.py" 40 "$A"; then
    $PY -u "$D/src/40_event_catalogue.py" "$A" \
        --nullfile "$D/results/permnull_40_$A.npz" || FAIL+=("40_$A:observed")
  else FAIL+=("40_$A:merge"); fi
  for MD in 4 6; do
    TAG=""; [[ $MD == 6 ]] && TAG="_d6"
    echo "################################################## 43 $A maxd=$MD"
    if $PY -u "$D/src/46_perm_merge.py" 43 "$A" --tag "$TAG"; then
      $PY -u "$D/src/43_soft_events.py" "$A" --maxd "$MD" \
          --nullfile "$D/results/permnull_43_${A}${TAG}.npz" || FAIL+=("43d${MD}_$A:observed")
    else FAIL+=("43d${MD}_$A:merge"); fi
  done
done
if [[ ${#FAIL[@]} -gt 0 ]]; then echo "FAILED: ${FAIL[*]}"; exit 1; fi
echo "ALL 15 CONFIGURATIONS COLLECTED"

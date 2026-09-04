#!/usr/bin/env bash
# Pool the low-floor permutation parts and run each observed scan against them.
# Not set -e: one arm failing must not stop the other fourteen.
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
LAM_HARD=2; LAM_SOFT=4
FAIL=(); cd "$ROOT"
for A in Mouse3 Mouse1 Mouse2 Initial Subclone; do
  echo "################################################## 40 $A lam=$LAM_HARD"
  if $PY -u "$D/src/46_perm_merge.py" 40 "$A" --tag "_lam$LAM_HARD"; then
    $PY -u "$D/src/40_event_catalogue.py" "$A" --lam $LAM_HARD \
        --nullfile "$D/results/permnull_40_${A}_lam$LAM_HARD.npz" || FAIL+=("40_$A")
  else FAIL+=("40_$A:merge"); fi
  for MD in 4 6; do
    TAG=""; [[ $MD == 6 ]] && TAG="_d6"
    echo "################################################## 43 $A d$MD lam=$LAM_SOFT"
    if $PY -u "$D/src/46_perm_merge.py" 43 "$A" --tag "${TAG}_lam$LAM_SOFT"; then
      $PY -u "$D/src/43_soft_events.py" "$A" --maxd "$MD" --lam $LAM_SOFT \
          --nullfile "$D/results/permnull_43_${A}${TAG}_lam$LAM_SOFT.npz" \
          || FAIL+=("43d${MD}_$A")
    else FAIL+=("43d${MD}_$A:merge"); fi
  done
done
if [[ ${#FAIL[@]} -gt 0 ]]; then echo "FAILED: ${FAIL[*]}"; exit 1; fi
echo "ALL 15 LOW-FLOOR CONFIGURATIONS COLLECTED"

#!/usr/bin/env bash
# Re-run the OBSERVED scans against the already-stored low-floor nulls, choosing
# each stratum's threshold by an ABSOLUTE expected-false budget (null mean <= 2
# per stratum, so <= ~12 expected false candidates per arm) instead of a 5% rate.
# No permutations are redone -- permnull_*_lam{2,4}.npz already hold B=1000.
ROOT=/data1/choij10/justin/pando
D="$ROOT/analyses/2026-08_park-compatibility"
PY=/data1/choij10/justin/envs/pando/bin/python
BUD=2
FAIL=(); cd "$ROOT"
for A in Mouse3 Mouse1 Mouse2 Initial Subclone; do
  echo "################################################## 40 $A budget=$BUD"
  $PY -u "$D/src/40_event_catalogue.py" "$A" --lam 2 --budget $BUD \
      --nullfile "$D/results/permnull_40_${A}_lam2.npz" || FAIL+=("40_$A")
  for MD in 4 6; do
    TAG=""; [[ $MD == 6 ]] && TAG="_d6"
    echo "################################################## 43 $A d$MD budget=$BUD"
    $PY -u "$D/src/43_soft_events.py" "$A" --maxd "$MD" --lam 4 --budget $BUD \
        --nullfile "$D/results/permnull_43_${A}${TAG}_lam4.npz" || FAIL+=("43d${MD}_$A")
  done
done
if [[ ${#FAIL[@]} -gt 0 ]]; then echo "FAILED: ${FAIL[*]}"; exit 1; fi
echo "ALL 15 REPORTED AT AN ABSOLUTE FALSE-EVENT BUDGET"

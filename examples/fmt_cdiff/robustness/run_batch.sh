#!/usr/bin/env bash
# Normal-priority batch, 4-way, specificity configs first, skips configs already done.
# Run detached: setsid nohup ./run_batch.sh </dev/null &>/dev/null &
cd /home/cumbof/muODE/examples/fmt_cdiff/robustness
PY=/home/cumbof/.conda/envs/muode-run/bin/python

# priority order: refs + leave-one-out (the headline) first, then the sensitivity sweeps
{ echo configs/ref_bile_on.json
  echo configs/ref_bile_off.json
  ls configs/loo_*.json
  ls configs/dil_*.json configs/dose_*.json configs/time_*.json
} > order.txt

run() {
  lbl=$(basename "$1" .json)
  if [ -s "results/$lbl.json" ]; then echo "SKIP $lbl" >> batch.log; return; fi
  "$PY" robust_one.py "$1" >> batch.log 2>&1 || echo "FAIL $1" >> batch.log
}
export -f run; export PY

cat order.txt | xargs -P 4 -I{} bash -c 'run "$@"' _ {}
echo BATCH_DONE >> batch.log

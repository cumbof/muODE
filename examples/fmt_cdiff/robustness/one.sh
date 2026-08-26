#!/usr/bin/env bash
# Single-config timed test, fully detachable. Writes timestamps + result to one.log.
cd /home/cumbof/muODE/examples/fmt_cdiff/robustness || exit 3
PY=/home/cumbof/.conda/envs/muode-run/bin/python
{
  echo "START $(date +%T)"
  echo "python: $($PY -c 'import sys;print(sys.version.split()[0])' 2>&1)"
  echo "import muode: $($PY -c 'import muode,cobra;print(\"ok\",cobra.__version__)' 2>&1)"
  t0=$SECONDS
  timeout 1200 "$PY" robust_one.py configs/ref_bile_off.json
  rc=$?
  echo "END $(date +%T) rc=$rc elapsed=$((SECONDS-t0))s"
} > one.log 2>&1

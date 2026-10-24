#!/bin/bash
# usage: run_stream.sh BIN THREADS CORELIST REPDIR dir...   (sequential, pinned)
BIN=$1; T=$2; CORES=$3; REP=$4; shift 4
for d in "$@"; do
  cd "$REP/$d" || exit 1
  OMP_NUM_THREADS=$T OMP_SCHEDULE=static OMP_PROC_BIND=true OMP_PLACES=cores taskset -c $CORES /usr/bin/time -p $BIN > stdout.log 2> time.log
  echo "$(date +%T) done $REP/$d $(grep -E '^(real|user)' time.log | tr '\n' ' ')"
done

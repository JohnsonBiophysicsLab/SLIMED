#!/bin/bash
# writes results_<rep>.txt in the local format: name real X user Y flips=N
for rep in rep1 rep2 t4; do
  : > results_$rep.txt
  for d in $rep/*/; do n=$(basename $d); [ -f $d/time.log ] || continue
    flips=$(awk -F, 'NR>1 && $8==1' $d/inputEdgeFlips.csv 2>/dev/null | wc -l | tr -d ' ')
    echo "$n $(grep -E '^(real|user)' $d/time.log | tr '\n' ' ') flips=$flips $(grep 'Edge flips' $d/stdout.log | tail -1 | sed 's/.*flips: //')" >> results_$rep.txt
  done
done

#!/bin/bash
#SBATCH --job-name=slimed_depth_sweep
#SBATCH --account=mjohn218
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --time=02:30:00
#SBATCH --output=%x_%j.out
set -u
source /etc/profile.d/modules.sh 2>/dev/null || true
module load GCC/12.3.0 >/dev/null 2>&1
export LD_LIBRARY_PATH=/data/apps/extern/easybuild/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH:-}
module list 2>&1 | grep -io "gcc[^ ]*" | tr '\n' ' '; echo
cd $SLURM_SUBMIT_DIR
BIN=$SLURM_SUBMIT_DIR/../SLIMED/build-omp/bin/membrane_dynamics
echo "node $(hostname)  $(date)"; lscpu | grep -E 'Model name|^CPU\(s\)|Thread|MHz'; uptime
ls -la $BIN; ldd $BIN | grep -E "stdc|gomp|gsl"
ps -u $USER -o pid,psr,pcpu,etime,comm | grep -v "ps\|grep\|bash\|slurm" || true
# pre-flight: 50 steps must run, else abort the job
rm -rf preflight && cp -r rep1/d3_n1000 preflight && cd preflight && sed -i "s/^maxIterations = .*/maxIterations = 50/" input.params
OMP_NUM_THREADS=1 /usr/bin/time -p $BIN > stdout.log 2> time.log; cd ..
if ! grep -q "ITERATION:49" preflight/stdout.log; then echo "PRE-FLIGHT FAILED"; cat preflight/time.log; exit 1; fi
echo "pre-flight ok: $(grep -E '^(real|user)' preflight/time.log | tr '\n' ' ')"
# 1-thread series: rep1 on cores 0-3, rep2 on cores 24-27 (second socket), 8 streams at once
pids=()
while read i list; do
  bash run_stream.sh $BIN 1 $i      $SLURM_SUBMIT_DIR/rep1 $list > stream_rep1_$i.log 2>&1 & pids+=($!)
  bash run_stream.sh $BIN 1 $((24+i)) $SLURM_SUBMIT_DIR/rep2 $list > stream_rep2_$i.log 2>&1 & pids+=($!)
done < streams.txt
wait "${pids[@]}"
echo "1-thread series done $(date)"; cat stream_rep*_*.log
# 4-thread wall-time series: two at a time, 4 cores each
bash run_stream.sh $BIN 4 0-3   $SLURM_SUBMIT_DIR/t4 default_n1000 default_n2000 frozen_n1000 > stream_t4_a.log 2>&1 &
bash run_stream.sh $BIN 4 24-27 $SLURM_SUBMIT_DIR/t4 legacy_n1000 legacy_n2000 frozen_n2000 > stream_t4_b.log 2>&1 &
wait
echo "4-thread series done $(date)"; cat stream_t4_*.log

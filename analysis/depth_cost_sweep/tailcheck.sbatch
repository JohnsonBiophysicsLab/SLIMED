#!/bin/bash
#SBATCH --job-name=slimed_tailcheck
#SBATCH --account=mjohn218
#SBATCH --partition=parallel
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=48
#SBATCH --exclusive
#SBATCH --time=00:40:00
#SBATCH --output=%x_%j.out
set -u
source /etc/profile.d/modules.sh 2>/dev/null || true
module load GCC/12.3.0 >/dev/null 2>&1
export LD_LIBRARY_PATH=/data/apps/extern/easybuild/GCCcore/12.3.0/lib64:${LD_LIBRARY_PATH:-}
cd $SLURM_SUBMIT_DIR
BIN=$SLURM_SUBMIT_DIR/../SLIMED/build-omp/bin/membrane_dynamics
echo "node $(hostname) $(date)"; uptime
i=0
for v in d6555_notail d5_tail default_rep3 d5_rep3; do
  bash run_stream.sh $BIN 1 $i $SLURM_SUBMIT_DIR ${v}_n1000 > stream_${v}_a.log 2>&1 &
  bash run_stream.sh $BIN 1 $((24+i)) $SLURM_SUBMIT_DIR ${v}_n2000 > stream_${v}_b.log 2>&1 &
  i=$((i+1))
done
wait; echo "done $(date)"; cat stream_*.log

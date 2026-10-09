#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:05:00
#SBATCH --nodes=1      ##adjust
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=boost_usr_prod

# to avoid perl warning
export LC_CTYPE=en_US.UTF-8
export LC_ALL=en_US.UTF-8
# load modules
module purge
module load nvhpc/23.1
module load openmpi/4.1.4--nvhpc--23.1-cuda-11.8

# the code uses paths relative to set_run (./sc_compiled, ./results, ./initial_fields):
# run from set_run, also when the job is submitted from the main folder (sbatch set_run/go.sh)
cd "${SLURM_SUBMIT_DIR:-.}"
if [ -d ./set_run/sc_compiled ]; then cd ./set_run; fi

#if using HPC-SDK, CUDA-aware already enabled):
mpirun -n NUMTASKS --map-by socket ./sc_compiled/flow36

# submit script with sbatch

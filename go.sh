#!/bin/bash
#SBATCH --account="IscrB_BMO"
#SBATCH --job-name="flo36gpu_test"
#SBATCH --time=00:05:00
#SBATCH --ntasks=NUMTASKS   ## total MPI tasks = NYCPU*NZCPU (set by compile.sh), nodes = ntasks/4
#SBATCH --ntasks-per-node=4
#SBATCH --gres=gpu:4   ###4 GPUs per node on 4 MPI tasks
#SBATCH --output=test.out
#SBATCH --error=test.err
#SBATCH --partition=boost_usr_prod
#SBATCH --qos=boost_qos_dbg   ## debug QOS: max 8 nodes, 30 min

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
# 4 consecutive ranks per node (one y-group with NYCPU=4: xz2yz/yz2xz stay on NVLink);
# binder_leo.sh sets the InfiniBand card of each GPU (UCX_NET_DEVICES=mlx5_<local rank>:1)
chmod +x ./binder_leo.sh
mpirun -n NUMTASKS --map-by ppr:4:node ./binder_leo.sh ./sc_compiled/flow36
# profiling (as in MHIT36): one report per rank, report_<rank>.nsys-rep
#mpirun -n NUMTASKS --map-by ppr:4:node nsys profile -t cuda,nvtx,mpi,openacc --nic-metrics=true -o report_%q{OMPI_COMM_WORLD_RANK} ./binder_leo.sh ./sc_compiled/flow36

# submit script with sbatch

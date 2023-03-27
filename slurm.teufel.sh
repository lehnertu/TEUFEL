#!/bin/bash
#SBATCH --job-name=TEUFEL    # name of the job
#SBATCH --partition=defq     # partition to be used (defq, gpu or intel)
#SBATCH --time=8:00:00       # walltime (up to 96 hours)
#SBATCH --nodes=2            # number of nodes
#SBATCH --ntasks=4           # number of MPI tasks
#SBATCH --ntasks-per-node=2  # number of MPI tasks per node
#SBATCH --cpus-per-task=20   # Number of CPU cores per task

module load gcc
module load openmpi
module load zlib
module load hdf5

cd $HOME/teufel              # path where data is located

# Fix distribution of tasks to nodes as workaround for bug in slurm
# Proposed by Henrik Schulz (FWCI)
export SLURM_TASKS_PER_NODE="$((SLURM_NTASKS / SLURM_NNODES))$( for ((i=2; i<=$SLURM_NNODES; i++)); \
do printf ",$((SLURM_NTASKS / SLURM_NNODES))"; done )"

# List of allocated nodes
scontrol show hostname $SLURM_JOB_NODELIST | paste -d, -s > hostsfile_job_$SLURM_JOBID.txt

export OMP_NUM_THREADS=20
mpiexec ./build/teufel examples/elbe-u300_200pC.xml

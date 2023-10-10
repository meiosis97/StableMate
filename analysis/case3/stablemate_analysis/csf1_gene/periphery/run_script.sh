#!/bin/bash
# Created by the University of Melbourne job script generator for SLURM
# Tue Aug 16 2022 14:43:14 GMT+1000 (澳大利亚东部标准时间)

# Partition for the job:
#SBATCH --partition=physical

# Multithreaded (SMP) job: must run on one node 
#SBATCH --nodes=21
#SBATCH --ntasks-per-node=1

# The name of the job:
#SBATCH --job-name="stablemate"

# The project ID which this job should run under:
#SBATCH --account="punim1662"

# Maximum number of tasks/CPU cores used by the job:
#SBATCH --cpus-per-task=1

# The amount of memory in megabytes per process in the job:
#SBATCH --mem=8G

# The maximum running time of the job in days-hours:mins:sec
#SBATCH --time=0-5:0:00

# stdo and stde
#SBATCH -o MPI.out # STDOUT 
#SBATCH -e MPI.err # STDERR

# check that the script is launched with sbatch
if [ "x$SLURM_JOB_ID" == "x" ]; then
   echo "You need to submit your job to the queuing system with sbatch"
   exit 1
fi

# Run the job from the directory where it was launched (default)

##DO NOT ADD/EDIT BEYOND THIS LINE##
##Job monitor command to list the resource usage
my-job-stats -a -n -s

echo "------------------------------------------------------------"
echo "SLURM JOB ID: $SLURM_JOBID"
echo "$SLURM_NTASKS tasks"
echo "$SLURM_NTASKS_PER_NODE tasks per node"
echo "Running on nodes: $SLURM_NODELIST"
echo "------------------------------------------------------------"

module load r/4.2.0
export R_LIBS_USER='~/libs/R_libs/'
mpirun -n 1 R --slave -f ./script.R
my-job-stats -a -n -s

touch DONE 

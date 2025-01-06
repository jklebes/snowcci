#!/bin/bash 

#SBATCH --job-name=snowcci

# launch ntasks jobs
./FSM2 < nlst_${SLURM_ARRAY_TASK_ID}

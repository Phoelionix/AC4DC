#!/bin/bash
#SBATCH --nodes=1
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=12
#SBATCH --time=1-23:59:00
#SBATCH --array=1-40
#SBATCH --output=slurmies/Z-%x.%j.out

# Load required modules
module purge
#bash load_modules.sh
source load_modules.sh

echo "Running $SLURM_ARRAY_TASK_COUNT simulations on $SLURM_CPUS_ON_NODE CPU cores"
echo "End time: $SLURM_JOB_END_TIME"

./ac4dc input/_batches/batch_SH2_Gd/SH2_Gd-${SLURM_ARRAY_TASK_ID}.mol
echo "Done!"
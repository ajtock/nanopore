#!/bin/bash

#! Account to use (users can be part of one or more accounts).
#! Currently we have only two: 'bioinf' and 'it-support'.
#SBATCH -A bioinf

#! Partition to run on
#SBATCH -p production

#! Number of nodes to allocate
#SBATCH --nodes=2

#! Email notification for job conditions (e.g., START, END,FAIL, HOLD or none)
#SBATCH --mail-type=END,FAIL

#! Do not re-queue job after system fail
#SBATCH --no-requeue

#! Per-array-task-ID log file (NOTE: enclosing directory must exist)
#SBATCH -o logs/among_read_variation_plotting_%a.log

#! Per-array-task-ID error file (NOTE: enclosing directory must exist)
#SBATCH -e logs/among_read_variation_plotting_error_%a.txt

#! Number of CPUs. Default: 1
#SBATCH -c 5

#! RAM. Default: 1G
#SBATCH --mem=1G

#! Time in HH:MM:SS. Default: 14 days (currently)
#SBATCH -t 00:01:00

#! Array task IDs
#SBATCH -a 2-3

#! Get the relevant line of the TSV parameters file
PARAMS=$(cat slurm_params/among_read_variation_plotting_slurm_params.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#! Get each parameter
SAMPLE=$(echo $PARAMS | cut -d ' ' -f 1)
REF_BASE=$(echo $PARAMS | cut -d ' ' -f 2)
GENOME_BIN_SIZE=$(echo $PARAMS | cut -d ' ' -f 3)
GENOME_STEP_SIZE=$(echo $PARAMS | cut -d ' ' -f 4)
CONTEXT=$(echo $PARAMS | cut -d ' ' -f 5)
NA_MAX=$(echo $PARAMS | cut -d ' ' -f 6)
CPU_PC=$(echo $PARAMS | cut -d ' ' -f 7)
CHR_NAME=$(echo $PARAMS | cut -d ' ' -f 8)

#! Output some informative messages
echo "Context: $CONTEXT"
echo "Slurm array task ID: $SLURM_ARRAY_TASK_ID"
echo "Number of CPUs used: $SLURM_CPUS_PER_TASK"
echo "This job is running on:"
hostname

#! Execute
./among_read_variation_plotting_slurm.R \
  $SAMPLE $REF_BASE \
  $GENOME_BIN_SIZE $GENOME_STEP_SIZE \
  $CONTEXT $NA_MAX \
  $CPU_PC $CHR_NAME

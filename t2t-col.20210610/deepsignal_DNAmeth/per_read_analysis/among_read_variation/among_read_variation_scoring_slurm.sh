#!/bin/bash

#! Account to use (users can be part of one or more accounts).
#! Currently we have only two: 'bioinf' and 'it-support'.
#SBATCH -A bioinf

#! Partition to run on
#SBATCH -p production

#! Email notification for job conditions (e.g., START, END,FAIL, HOLD or none)
#SBATCH --mail-type=END,FAIL
#SBATCH --mail-user=ajt200@cam.ac.uk # Doesn't seem to send these to my email inbox

#! Do not re-queue job after system fail
#SBATCH --no-requeue

#! Per-array-task-ID log file (NOTE: enclosing directory must exist)
#SBATCH -o logs/among_read_variation_scoring_%a.log

#! Per-array-task-ID error file (NOTE: enclosing directory must exist)
#SBATCH -e logs/among_read_variation_scoring_error_%a.txt

#! Number of nodes to allocate
#SBATCH --nodes=3

#! Number of CPUs per task. Default: 1
#SBATCH --cpus-per-task=32

#! Minimum RAM needed for all tasks
#! NOTE: Doesn't work with > 1M ; e.g., with 10M, 100M, 1G, get these errors:
#! sbatch: error: Memory specification can not be satisfied
#! sbatch: error: Batch job submission failed: Requested node configuration is not available
#SBATCH --mem=127G

#! Time in HH:MM:SS. Default: 14 days (currently)
#SBATCH -t 99:00:00

#! Array task IDs
#SBATCH -a 2-4

#! Get the relevant line of the TSV parameters file
PARAMS=$(cat slurm_params/among_read_variation_scoring_slurm_params.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

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
./among_read_variation_scoring.R \
  $SAMPLE $REF_BASE \
  $GENOME_BIN_SIZE $GENOME_STEP_SIZE \
  $CONTEXT $NA_MAX \
  $CPU_PC $CHR_NAME

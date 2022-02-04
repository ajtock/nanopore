#!/bin/bash

#SBATCH -D /home/ajt200/analysis/nanopore/t2t-col.20210610/deepsignal_DNAmeth

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
#SBATCH -o logs/per_read_tsv2bed_bam_conversion_per_chr_%a.log

#! Per-array-task-ID error file (NOTE: enclosing directory must exist)
#SBATCH -e logs/per_read_tsv2bed_bam_conversion_per_chr_error_%a.txt

#! Number of nodes to allocate
#SBATCH --nodes=1

#! Number of CPUs per task. Default: 1
#SBATCH --cpus-per-task=32

#! Minimum RAM needed for all tasks
#! NOTE: Doesn't work with > 1M ; e.g., with 10M, 100M, 1G, get these errors:
#! sbatch: error: Memory specification can not be satisfied
#! sbatch: error: Batch job submission failed: Requested node configuration is not available
#SBATCH --mem=128G

#! Time in HH:MM:SS. Default: 14 days (currently)
#SBATCH -t 336:00:00

#! Array task IDs
#SBATCH -a 2-3

#! Get the relevant line of the TSV parameters file
PARAMS=$(cat slurm_params/per_read_tsv2bed_bam_conversion_per_chr_slurm_params.tsv | head -n $SLURM_ARRAY_TASK_ID | tail -n 1)

#! Get each parameter
PREFIX=$(echo $PARAMS | cut -d ' ' -f 1)
REF_BASE=$(echo $PARAMS | cut -d ' ' -f 2)
THREADS=$(echo $PARAMS | cut -d ' ' -f 3)
SORT_MEM=$(echo $PARAMS | cut -d ' ' -f 4)
CONTEXT=$(echo $PARAMS | cut -d ' ' -f 5)

#! Output some informative messages
echo "Context: $CONTEXT"
echo "Slurm array task ID: $SLURM_ARRAY_TASK_ID"
echo "Number of CPUs used: $SLURM_CPUS_PER_TASK"
echo "This job is running on:"
hostname
echo $(which tabix)
echo $(which samtools)

#! Execute
./per_read_tsv2bed_bam_conversion_per_chr.sh \
  $PREFIX \
  $REF_BASE \
  $THREADS \
  $SORT_MEM \
  $CONTEXT

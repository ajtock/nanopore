#!/bin/bash

#! Account to use (users can be part of one or more accounts).
#! Currently we have only two: 'bioinf' and 'it-support'.
#SBATCH -A bioinf

#! Partition to run on
#SBATCH -p production

#! Number of nodes to allocate
#SBATCH --nodes=1

#! Number of CPUs per task. Default: 1
#SBATCH --cpus-per-task=2

#! Minimum RAM needed for all tasks
#SBATCH --mem=1G

#! Time in HH:MM:SS. Default: 14 days (currently)
#SBATCH -t 00:01:00

#! log file (enclosing directory must exist)
#SBATCH -o logs/simple_job.log

#! error file (enclosing directory must exist)
#SBATCH -e logs/simple_job_error.txt


source activate ChIPseq_mapping
which samtools

#! Output some informative messages
echo "Number of CPUs used: $SLURM_CPUS_PER_TASK"
echo "This job is running on:"
hostname

conda deactivate

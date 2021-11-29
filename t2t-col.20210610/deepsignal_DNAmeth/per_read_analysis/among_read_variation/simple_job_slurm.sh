#!/bin/bash

#! Working directory
#!SBATCH -D /home/ajt200/analysis/nanopore/t2t-col.20210610/deepsignal_DNAmeth/per_read_analysis/among_read_variation

#! Account to use (users can be part of one or more accounts).
#! Currently we have only two: 'bioinf' and 'it-support'.
#SBATCH -A bioinf

#! Per-array-task-ID log file (enclosing directory must exist)
#SBATCH -o logs/simple_job.log

#! Per-array-task-ID error file (enclosing directory must exist)
#SBATCH -e logs/simple_job_error.txt

#! Output some informative messages
echo "Number of CPUs used: $SLURM_CPUS_PER_TASK"
echo "This job is running on:"
hostname

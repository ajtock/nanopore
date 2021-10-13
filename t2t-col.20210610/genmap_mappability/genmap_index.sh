#!/bin/bash

# Build a reference genome index for mappability computation
# using GenMap version 1.3.0
# "GenMap is a tool for fast and exact computation of genome mappability
# and can also be used for multiple genomes, e.g., to search for marker sequences."

# Usage on node7:
# /scripts/csmit -m 20G -c 1 "bash genmap_index.sh t2t-col.20210610.fa"

genome=$1
toolDir=/home/ajt200/miniconda3/envs/ChIPseq_mapping/bin

${toolDir}/genmap index --fasta-file ${genome} \
                        --index index/ \
                        --algorithm divsufsort

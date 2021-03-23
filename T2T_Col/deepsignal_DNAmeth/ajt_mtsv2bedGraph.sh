#!/bin/bash

# Usage on node7 (must be allocated to node10 as requires > 300 GB RAM):
# csmit -m 400G -c 1 "bash ./ajt_mtsv2bedGraph.sh WT_deepsignalDNAmeth_95_30kb_MappedOn_T2T_Col CpG" 

prefix=$1
context=$2

source activate ChIPseq_mapping

ajt_mtsv2bedGraph.py --input ${prefix}_${context}_raw.tsv \
                     --mod ${context} \
                     --call-threshold 0.02 \
                     --genome /home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa \
                     --offset 0 \
                     --verbose | \
  sort -k1,1 -k2,2n | bgzip > ${prefix}_${context}.bed.gz
tabix -p bed ${prefix}_${context}.bed.gz

conda deactivate

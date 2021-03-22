#!/bin/bash

# Usage on node7:
# csmit -m 20G -c 1 "bash ./ajt_mtsv2bedGraph.sh head_WT_deepsignalDNAmeth_95_30kb_MappedOn_T2T_Col CpG" 

prefix=$1
context=$2

#/home/nm359/nanopore-methylation-utilities/mtsv2bedGraph.py -i ${filename}.tsv -c ${threshold} -g ${reference} -v |\
#  sort -k1,1 -k2,2n | bgzip > ${output}.bed.gz
#tabix -p bed ${output}.bed.gz

#source activate ChIPseq_mapping

ajt_mtsv2bedGraph.py --input ${prefix}_${context}_raw.tsv \
                     --mod ${context} \
                     --call-threshold 0.02 \
                     --genome /home/ajt200/analysis/nanopore/T2T_Col/T2T_Col.fa \
                     --offset 0 \
                     --verbose > ${prefix}_${context}.bed
#  sort -k1,1 -k2,2n | bgzip > ${prefix}.bed.gz
#tabix -p bed ${prefix}.bed.gz

#conda deactivate

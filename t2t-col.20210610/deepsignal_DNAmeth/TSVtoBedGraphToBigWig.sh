#!/bin/bash

# Convert TSV file containing deepsignal-calculated
# context-specific DNA methylation frequencies into bigWig format
# suitable for use with deepTools computeMatrix for fine-scale profiling

# Usage:
# ./TSVtoBedGraphToBigWig.sh Col_0_deepsignalDNAmeth_30kb_90pc_MappedOn_t2t-col.20210610_CpG

i=$1

source activate BSseq_mapping

[ -d ./bg ] || mkdir ./bg
[ -d ./bw ] || mkdir ./bw

# USAGE: bedGraphToBigWig in.bedGraph chrom.sizes out.bw                                                  
# where in.bedGraph is a four-column file in the format:                                                  
#       <chrom> <start> <end> <value>
# and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>                             
# and out.bw is the output indexed big wig file.
# The input bedGraph file must be sorted, use the unix sort command:                                      
(cat ${i}.tsv | LC_COLLATE=C sort -k1,1 -k2,2n \
| awk 'BEGIN{FS="\t";OFS="\t"} {print $1, $2, $2+1, $10*100}' - \
> bg/${i}.bedGraph; \
bedGraphToBigWig bg/${i}.bedGraph /home/ajt200/analysis/nanopore/t2t-col.20210610/t2t-col.20210610.fa.sizes bw/${i}.bw ) \
&> ${i}_TSVtoBedGraphToBigWig.log

source deactivate

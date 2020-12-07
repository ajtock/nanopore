#!/bin/bash

# Compute mappability for genmap-indexed reference genome (see genmap_index.sh) 
# using GenMap version 1.3.0
# "GenMap is a tool for fast and exact computation of genome mappability
# and can also be used for multiple genomes, e.g., to search for marker sequences."

# Usage on node7:
# /scripts/csmit -m 200G -c 1 "bash genmap_map.sh 40 2 Athaliana_ONT_RaGOO_v2.0"

K=$1
E=$2
prefix=$3

[ -d "map_K"${K}"_E"${E} ] || mkdir "map_K"${K}"_E"${E}

source activate ChIPseq_mapping

genmap map --length ${K} \
           --errors ${E} \
           --index index \
           --output "map_K"${K}"_E"${E} \
           --txt --wig --bedgraph

conda deactivate

# Convert genmap-generated bedGraph file into bigWig format
# suitable for use with deepTools computeMatrix for fine-scale profiling

source activate BSseq_mapping

# USAGE: bedGraphToBigWig in.bedGraph chrom.sizes out.bw                                                  
# where in.bedGraph is a four-column file in the format:                                                  
#       <chrom> <start> <end> <value>
# and chrom.sizes is a two-column file/URL: <chromosome name> <size in bases>                             
# and out.bw is the output indexed big wig file.
# The input bedGraph file must be sorted; use the unix sort command:                                      
(cat "map_K"${K}"_E"${E}/${prefix}.genmap.bedgraph \
| LC_COLLATE=C sort -k1,1 -k2,2n > "map_K"${K}"_E"${E}/${prefix}.genmap.sorted.bedgraph
bedGraphToBigWig "map_K"${K}"_E"${E}/${prefix}.genmap.sorted.bedgraph \
                 "map_K"${K}"_E"${E}/${prefix}.genmap.chrom.sizes \
                 "map_K"${K}"_E"${E}/${prefix}.genmap.bw ) \
&> "map_K"${K}"_E"${E}/${prefix}_bedGraphToBigWig.log

conda deactivate

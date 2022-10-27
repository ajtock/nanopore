#!/bin/bash


source activate R-4.0.0
./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_H3K9me2_Rep1_ChIP' 'WT_H3K9me2_Rep1_input' '170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610' '170101_Chris_H3K9me2_ChIP/WT/snakemake_ChIPseq_t2t-col.20210610' 'H3K9me2' 'H3K9me2 input' 'dodgerblue1' 'dodgerblue1' 't2t-col.20210610'
./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_REC8_HA_Rep2_ChIP' 'WT_REC8_Myc_Rep1_input' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8' 'REC8 input' 'firebrick1' 'firebrick1' 't2t-col.20210610'
./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_ASY1_Rep1_ChIP' 'WT_REC8_Myc_Rep1_input' '20190722_cal66_Athaliana_ChIPseq_ASY1/fastq_pooled/snakemake_ChIPseq_t2t-col.20210610' 'REC8_pooled/snakemake_ChIPseq_t2t-col.20210610' 'ASY1' 'REC8 input' 'darkgreen' 'firebrick1' 't2t-col.20210610'
./LRZgenes_NRZgenes_armgenes_ChIP_1metaprofile.R 'Chr1,Chr2,Chr3,Chr4,Chr5' both 2000 2000 '2kb' 10 10bp '0.02,0.96' 'WT_SPO11oligos_Rep1' 'WT_gDNA_Rep1_R1' '160518_Kyuha_SPO11oligos/WT/snakemake_SPO11oligos_t2t-col.20210610' '150701_Natasha_gDNA/WT/R1/snakemake_SPO11oligos_t2t-col.20210610' 'SPO11-1-oligos' 'gDNA' 'magenta3' 'magenta3' 't2t-col.20210610'
conda deactivate



#!/bin/bash

source activate ChIPseq_mapping
snakemake -p --cores 48
conda deactivate

#!/bin/bash

# Usage:
# ./convert_gene_Liftoff_GFF3_to_prot_fa.sh t2t-col.20210610

genome=$1

/applications/cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread ${genome}.genes.gff3 \
                                                             -g ${genome}.fa \
                                                             -y ${genome}.genes.prot.fa

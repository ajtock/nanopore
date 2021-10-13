#!/bin/bash

# Usage:
# ./convert_EDTA_Chr3_GFF3_to_GTF.sh t2t-col.20210610

genome=$1

/applications/cufflinks/cufflinks-2.2.1.Linux_x86_64/gffread ${genome}.fasta.mod.EDTA.TEanno.Chr3.gff3 \
                                                             -T -o ${genome}.fasta.mod.EDTA.TEanno.Chr3.gtf

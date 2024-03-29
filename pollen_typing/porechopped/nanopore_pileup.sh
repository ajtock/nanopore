#!/bin/bash

# Using samtools (version 1.9) mpileup, create a variant pileup file
# in which each column of variant alleles corresponds to those
# called within one ONT read alignment

# The resultant pileup file can then be analysed using pileup_to_haplotype.R to
# evaluate recombination rate (cM/Mb) at each marker

# Usage:
# bash ./nanopore_pileup.sh barcode1 60 '3:634109-639934' 48

#sample=barcode1
#MAPQ=60
#region='3:634109-639934'
#threads=48
 
sample=$1
MAPQ=$2
region=$3
threads=$4

# Activate conda software enrivonment
source activate longshot
which samtools 
#/home/ajt200/anaconda3/envs/longshot/bin/samtools
samtools --version
#samtools 1.9
#Using htslib 1.9

# Create BED3 file containing known variant sites of interest
# to be proovided to mpileup.
# Need to create this once only (e.g., for first sample)
if [ ${sample} == "barcode1" ]
then
  awk 'BEGIN {OFS="\t"}; {print $1, ($2 - 1), $2}' 3a_polymorphisms/3a_SNPs_indels_in_Col_and_Ws.tsv | \
    tail -n +2 > 3a_polymorphisms/3a_SNPs_indels_in_Col_and_Ws.bed
fi

# Create per-alignment BAM output directory if it doesn't already exist
[[ -d aln ]] || mkdir aln
[[ -d aln/${sample} ]] || mkdir aln/${sample}

# Filter BAM to retain only primary alignments and
# exclude alignments with MAPQ scores < ${MAPQ}
# on Minimap2-assigned MAPQ of 60, see https://github.com/lh3/minimap2/issues/447
# and https://github.com/lh3/minimap2/issues/528
# -F 2304 excludes secondary and supplementary alignments
samtools view -bh -@ ${threads} -F 2304 -q ${MAPQ} \
  -o bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam \
  bam/MM2_TAIR10_${sample}_sorted.bam
samtools index bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam 

# Extract SAM header for later concatenation with each alignment in the BAM file.
# The filenames of these individual alignments in BAM format will be listed
# in a TXT file to be supplied to samtools mpileup
samtools view -H bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam \
  > bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_header.sam

# Get the number of alignments in the BAM file
aln=$(samtools view bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam | wc -l)
echo ${aln}
#78302

# Generate indexed BAM file for each alignment
# Method 1 (much more efficient than Method 2)
# IFS="" preserves leading and trailing white space in $LINE,
# and prevents the read function from splitting up each line into fields
# The -r option disables the interpretation of backslashes as escape sequences
# Because the read function fails when it encounters end-of-file before the line ends,
# read -r line || [[ -n "$LINE" ]] tests for a non-empty line
# The 4th line of this while loop was previously divided between
# multiple lines (using " \" after each "|")
# This works when run in an interactive bash shell but for
# some reason fails when run as a bash script
counter=0
while IFS="" read -r line || [[ -n "${line}" ]]
do
  ((counter++))
  printf '%s\n' "${line}" | /bin/cat bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_header.sam - | samtools view -bh > aln/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_${counter}.bam
done < <(samtools view bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam)

# Sanity check to ensure that the number of alignments in the original BAM file
# is the same as the number of separate per-alignment BAM files in aln/${sample}
[[ $(( ${aln} )) == $( ls aln/${sample}/ | wc -l ) ]] || { echo >&2 "number of alignments in bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam is NOT equal to the number of separate per-alignment BAM files in aln/${sample}/ !"; exit 1; }

## Generate indexed BAM file for each alignment
## Method 2
## ${x} in this for loop decreases from ${aln} to 1
#for x in $(seq ${aln} -1 1)
#do
#  samtools view MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam | \
#    head -n $(( -${x} + 1)) | tail -n 1 | \
#    /bin/cat MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_header.sam - | \
#    samtools view -bh > aln/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_"$(( ${aln} - ( ${x} - 1 ) ))".bam
#  samtools index aln/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_"$(( ${aln} - ( ${x} - 1 ) ))".bam
#done

# Extract the total number of mismatches and gaps (NM:i:# SAM field)
# and the number of ambiguous bases in each alignment (nn:i:# SAM field)
counter=0
while IFS="" read -r line || [[ -n "$LINE" ]]
do
  ((counter++))
  printf '%s\n' "${line}" | awk 'BEGIN {OFS="\t"}; {print $12, $15}' - >> bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_NM_nn.txt
done < <(samtools view bam/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.bam)

# Generate list of per-alignment BAM file names (as generated by above loop)
# as a text file with one BAM file name per line
[[ -f aln/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_list.txt ]] && \
echo "ERROR: aln/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_list.txt EXISTS; NOT OVERWRITTEN" || \
for x in $(seq 1 ${aln})
do
  echo "aln/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_${x}.bam" \
    >> aln/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_list.txt
done

# Change limit on number of open files from 1024 to 1048576
ulimit -Sn 1048576

# Create variant pileup directory if it doesn't already exist
[[ -d plp ]] || mkdir plp

# Generate textual pileup matrix with columns (info in columns 4-6
# is output for each BAM file supplied via the --bam_list option):
# 1. Sequence name (e.g., chromosome)
# 2. 1-based coordinate
# 3. Reference genome base
# 4. BAM file 1: Number of reads covering this position
# 5. BAM file 1: Read bases at this position
## ("." denotes a match to the reference base on the positive strand,
##  "," denotes a match to the reference base on the negative strand,
##  uppercase letters denote non-reference read base on positive strand,
##  lowercase letters denote non-reference read base on negative strand,
##  "*" denotes the absence of a base call (e.g., due to no coverage),
##  "^]" denotes a base at the start of a read,
##  "$" denotes a base at the end of a  read)
# 6. BAM file 1: Phred-scaled base qualities
## see https://davetang.org/muse/2015/08/26/samtools-mpileup/
## http://www.htslib.org/doc/1.9/samtools.html
[[ -f plp/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.plp ]] && \
echo "ERROR: plp/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.plp EXISTS; NOT OVERWRITTEN" || \
samtools mpileup --bam-list aln/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_list.txt \
                 --fasta-ref TAIR10_chr_all.fa \
                 --positions 3a_polymorphisms/3a_SNPs_indels_in_Col_and_Ws.bed \
		 --ignore-RG \
                 --min-BQ 0 \
                 > plp/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}.plp

## Create per-alignment VCF output directory if it doesn't already exist
#[[ -d vcf ]] || mkdir vcf
#[[ -d vcf/${sample} ]] || mkdir vcf/${sample}
## Provide each aligment in BAM as input to longshot
#for x in $(seq 1 ${aln})
#do
#  longshot --region ${region} \
#           --bam aln/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_${x}.bam \
#           --ref TAIR10_chr_all.fa \
#           --out vcf/${sample}/MM2_TAIR10_${sample}_sorted_primary_MAPQ${MAPQ}_aln_${x}.vcf \
#           --min_cov 1 \
#           --min_alt_count 1 \
#           --min_allele_qual 1.0 \
#           --potential_snv_cutoff 1000 \
#           --hap_converge_delta 0.1
#done

source deactivate

#!/bin/bash

#This will extract all ONT read IDs from a gzip-compressed FASTQ file, assuming that the ONT read header lines each contain the pattern "read="

grep "read=" reads.fastq | sed 's/ runid.\+//g' - > readIDs.txt
zcat reads.fastq.gz | grep "read=" | sed 's/ runid.\+//g' - > readIDs.txt

#Once you have your read IDs of interest you could use an existing tool like seqtk subseq or a grep command that iterates over the lines in your readIDs.txt file:

grep -A 3 -f readIDs.txt reads.fastq > reads_subset.fastq

# Fastest
grep -A 3 -f readIDs.txt barcode4_2.fastq > reads_subset.fastq
zcat barcode4_2.fastq.gz | grep -A 3 -f readIDs.txt > reads_subset.fastq

# IFS="" preserves leading and trailing white space in $LINE,
# and prevents the read function from splitting up each line into fields
# The -r option disables the interpretation of backslashes as escape sequences
# Because the read function fails when it encounters end-of-file before the line ends,
# read -r line || [[ -n "$LINE" ]] tests for a non-empty line
cat readIDs.txt | while IFS="" read -r line || [[ -n "${line}" ]]
do
  zcat reads.fastq.gz | grep -A 3 "${line}"
done > reads_subset.fastq
gzip --best reads_subset.fastq

# Or a simpler for loop without the above constraints imposed:
for line in $(cat readIDs.txt)
do
  zcat reads.fastq.gz | grep -A 3 "${line}"
done > reads_subset.fastq
gzip --best reads_subset.fastq

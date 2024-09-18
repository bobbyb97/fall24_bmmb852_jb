#!/bin/bash

# Defining variables
GENOME="your_accession_here" #Example: GENOME=GCF_000188095.3

##--NOTHING BELOW THIS LINE SHOULD BE MODIFIED--##

if [ $GENOME == "your_accession_here" ]; then
    echo "ERROR: Make sure you define the 'GENOME' variable."
    exit 1
fi

# Test to see if bioinfo env is loaded
if command -v datasets; then 
# Download gff3 file
datasets download genome accession ${GENOME} --include genome,gff3,cds

# unzip the file
unzip ncbi_dataset.zip

# Separate genes into a different file.
cat ncbi_dataset/data/${GENOME}/*.gff | awk '$3 == "gene"' > ${GENOME}_genes.gff

# Quick count of annotated genes
echo "Number of annotated genes:"
cat ${GENOME}_genes.gff|wc -l

# Summary of most frequent annotation features
echo "Top 10 most frequent annotation features:"
cat ncbi_dataset/data/${GENOME}/*.gff | grep -v "^#" | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 10

else
    echo "ERROR: Make sure you run 'conda activate bioinfo' before running this script."
    exit 1
fi
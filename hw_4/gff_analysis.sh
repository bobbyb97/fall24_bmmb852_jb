#!/bin/bash

# This script will download a gff file from NCBI and print out a basic summary


# Defining variables
GENOME="GCF_000188095.3" # Replace with your desired accession
OUT_DIR=~ # Default output directory is home directory, change if desired

##--NOTHING BELOW THIS LINE SHOULD BE MODIFIED--##

if [ $GENOME == "your_accession_here" ]; then
    echo "ERROR: Make sure you define the variables"
    exit 1
fi

# Test to see if bioinfo environment is loaded
if command -v datasets; then 

# Create a bmmb_852_tmp directory to store files
echo "Creating a temporary directory to store files..."
mkdir ${OUT_DIR}/bmmb_852_tmp

# Download gff3 file
datasets download genome accession ${GENOME} --include genome,gff3,cds --filename ${OUT_DIR}/bmmb_852_tmp/ncbi_dataset.zip

# unzip the file
unzip ${OUT_DIR}/bmmb_852_tmp/ncbi_dataset.zip -d ${OUT_DIR}/bmmb_852_tmp

# Separate genes into a different file.
cat ${OUT_DIR}/bmmb_852_tmp/ncbi_dataset/data/${GENOME}/*.gff | awk '$3 == "gene"' > ${OUT_DIR}/bmmb_852_tmp/${GENOME}_genes.gff

# Quick count of annotated genes
echo "Number of annotated genes:"
cat ${OUT_DIR}/bmmb_852_tmp/${GENOME}_genes.gff|wc -l

# Summary of most frequent annotation features
echo "Top 10 most frequent annotation features:"
cat ${OUT_DIR}/bmmb_852_tmp/ncbi_dataset/data/${GENOME}/*.gff | grep -v "^#" | awk '{print $3}' | sort | uniq -c | sort -nr | head -n 10

# Remove the temporary directory
echo "Removing the temporary directory...done"
rm -rf ${OUT_DIR}/bmmb_852_tmp

else
    echo "ERROR: Make sure you run 'conda activate bioinfo' before running this script."
    exit 1
fi
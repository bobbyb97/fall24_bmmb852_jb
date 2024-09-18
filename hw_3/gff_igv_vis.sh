#!/bin/bash

# Download gff3 file
datasets download genome accession ${1} --include genome,gff3

# unzip the file
unzip ncbi_dataset.zip

# Separate intervals of type "gene" into a different file.

cat ncbi_dataset/data/${1}/*.gff | awk '$3 == "gene"' > ncbi_dataset/data/${1}/genes.gff

# Quick count of annotated genes
echo "Number of annotated genes:"
cat ncbi_dataset/data/${1}/genes.gff|wc -l


#Github link
#https://github.com/bobbyb97/fall24_bmmb852_jb/blob/main/hw_3/gff_igv_vis.md
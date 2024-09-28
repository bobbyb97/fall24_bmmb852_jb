#!bin/bash

set -uex

# Accession
ACC="GCA_036780745.1"

# Species
SPECIES="V_bombi"

# Number of reads
N=240000

# Length of reads
L=100

# Read files
R1=reads/wgsim_read1.fq
R2=reads/wgsim_read2.fq

##--NOTHING SHOULD BE MODIFIED BELOW THIS LINE--##

# Download the genome
# First check to see if it already exists

if [ -f ncbi_dataset/data/${ACC}/*.fna ]; then
    echo "Genome already downloaded"
else
    echo "Downloading genome"
    # Command to download genome
    datasets download genome accession ${ACC}
    # Unzip files, overwrite if necessary
    unzip -o ncbi_dataset.zip
    # Rename file
    ln -sf ncbi_dataset/data/${ACC}/*.fna ${SPECIES}.fna
fi

## Generate summary stats for your genome ##

# Size of file
wc -c ${SPECIES}.fna

# Size of genome and number of sequences
seqkit stats ${SPECIES}.fna

# Number of contigs
echo "Number of contigs:"
seqkit seq -i ${SPECIES}.fna | grep ">" | wc -l

# Get each contig ID and length
echo "Summary of contigs:"
seqkit fx2tab -n -l ${SPECIES}.fna | awk '{print $1 "\t" $NF}'


## Generate a simulated fastq file with wgsim ##

# Make the directory that will hold the reads extracts 
# the directory portion of the file path from the read
mkdir -p reads

# Simulate with no errors and no mutations
wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${SPECIES}.fna ${R1} ${R2}

# Run read statistics
seqkit stats ${R1} ${R2}
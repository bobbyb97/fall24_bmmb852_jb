#!/bin/bash

# Define variables
SRA="ERR12869977"
output_dir="reads"

##--FASTQ DOWNLOAD AND QC--##

# Check if bioinfo environment is loaded
if [ -z "$CONDA_DEFAULT_ENV" ]; then
    echo "Please load bioinfo environment"
    exit 1
fi

# Define read names
read1="${SRA}_1.fastq"
read2="${SRA}_2.fastq"

# Download first 10000 reads using fastq-dump
echo "Downloading from ${SRA} to ${output_dir} directory"

fastq-dump -X 10000 -F --outdir reads --split-files ${SRA}

# QC using fastqc 
echo "Running fastqc on ${output_dir}/${read1} and ${output_dir}/${read2}"
mkdir -p ${output_dir}/fastqc
fastqc ${output_dir}/${read1} ${output_dir}/${read2} -o ${output_dir}/fastqc

# Trim adapters and low quality bases using fastp
echo "Trimming adapters and low quality bases"
mkdir -p ${output_dir}/trimmed_fastq

  fastp --thread 4 \
        --in1 "${output_dir}/${read1}" \
        --in2 "${output_dir}/${read2}" \
        --out1 "${output_dir}/trimmed_fastq/${read1}_trimmed.fq.gz" \
        --out2 "${output_dir}/trimmed_fastq/${read2}_trimmed.fq.gz" \
        -g \
       	-q 35 \
       	-u 20 

# QC again
echo "Running fastqc on trimmed reads"
mkdir -p ${output_dir}/trimmed_fastq/fastqc

fastqc ${output_dir}/trimmed_fastq/${read1}_trimmed.fq.gz \
       ${output_dir}/trimmed_fastq/${read2}_trimmed.fq.gz \
       -o ${output_dir}/trimmed_fastq/fastqc 
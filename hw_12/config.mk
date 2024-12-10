### Variables ###
ACC ?= GCA_000287275.1
SPECIES ?= Carsonella_ruddii
SRA ?= DRR591238


##-- NOTHING BELOW THIS LINE SHOULD BE MODIFIED --##
REF = ${DOWNLOAD_DIR}/${SPECIES}.fna
GFF = ${DOWNLOAD_DIR}/${SPECIES}.gff
SCRIPT_DIR = make_vcf

## download.mk ##
DOWNLOAD_DIR = ${SPECIES}/data
SRA_DIR = ${DOWNLOAD_DIR}/${SRA}_reads
read1 = ${SRA}_1.fastq.gz
read2 = ${SRA}_2.fastq.gz

##fastp.mk##
TRIM_DIR = ${DOWNLOAD_DIR}/trimmed_fastq
trimmed_read1 = ${SRA}_trimmed_1.fq.gz
trimmed_read2 = ${SRA}_trimmed_2.fq.gz

##index.mk##
SORTED_BAM = ${BAM_DIR}/${SPECIES}_${SRA}_sorted.bam

##align.mk##

##bwa.mk##
BAM_DIR = ${DOWNLOAD_DIR}/bam
BAM = ${BAM_DIR}/${SPECIES}_${SRA}.bam
FLAGSTAT = ${BAM}_flagstat.txt

##vcf.mk##
VCF_DIR = ${DOWNLOAD_DIR}/vcf
IDX = ${DOWNLOAD_DIR}/idx/snpEff
LABEL = genome
SNPEFF_DB ?= ${IDX}/${LABEL}/snpEffectPredictor.bin

HTML = results/snpeff.html
CSV = results/snpeff.html
VCF = ${DOWNLOAD_DIR}/vcf/${SPECIES}.vcf
EFF = results/snpeff.vcf.gz


### Make setup ###

# Set the shell the commands run in.
SHELL := bash

# Execute all commands in a single shell.
.ONESHELL:

# Delete target files if the command fails.
.DELETE_ON_ERROR:

# Run the shell with strict error checking.
.SHELLFLAGS := -eu -o pipefail -c

# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables

# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules

# Define phony targets.
.PHONY: usage clean run
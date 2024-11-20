## Variables ##
ACC = GCA_000287275.1
SPECIES = Carsonella_ruddii
SRA = DRR591238

N = 1000
L = 100
R1 = reads/wgsim_read1.fq
R2 = reads/wgsim_read2.fq

output_dir = reads
read1 = $(SRA)_1.fastq
read2 = $(SRA)_2.fastq
trimmed_read1 = $(read1)_trimmed.fq.gz
trimmed_read2 = $(read2)_trimmed.fq.gz

IDX = idx/snpEff
LABEL = genome
SNPEFF_DB ?= ${IDX}/${LABEL}/snpEffectPredictor.bin
REF = $(SPECIES).fna
GFF = $(SPECIES).gff

HTML = results/snpeff.html
CSV = results/snpeff.html
VCF = vcf/${SPECIES}.vcf
EFF = results/snpeff.vcf.gz

# Set the shell the commands run in.
SHELL = bash

# Execute all commands in a single shell.
.ONESHELL:

# Run the shell with strict error checking.
.SHELLFLAGS = -eu -o pipefail -c

# Delete target files if the command fails.
.DELETE_ON_ERROR:

# Warn if a variable is not defined.
MAKEFLAGS += --warn-undefined-variables

# Disable built-in rules.
MAKEFLAGS += --no-builtin-rules

# Define phony targets.
.PHONY: usage conda_env genome simulate download trim clean all

## GENERAL ##
usage:
	@echo "usage: make [target] [OPTIONS]"
	@echo ""
	@echo "*** TARGETS ***"
	@echo ""
	@echo "General:"
	@echo "---------------"
	@echo "  usage		- Print the help"
	@echo "  conda_env	- Check if bioinfo environment is loaded"
	@echo "  variables	- Print optional variables, will automatically print defaults"
	@echo "  all		- Run all targets"
	@echo "  clean		- Remove all files"
	@echo ""
	@echo "Download:"
	@echo "---------------"
	@echo "  genome	- Download the genome"
	@echo "  download	- Download reads from SRA"
	@echo ""
	@echo "QC:"
	@echo "---------------"
	@echo "  trim		- Trim reads"
	@echo "  flagstat	- Run samtools flagstat"
	@echo "  filter	- Filter reads"
	@echo ""
	@echo "Align:"
	@echo "---------------"
	@echo "  index		- Index the genome"
	@echo "  align		- Align reads to the genome"
	@echo ""
	@echo "Simulate:"
	@echo "---------------"
	@echo "  simulate	- Simulate reads for the genome"
	@echo "  simalign	- Align simulated reads to the genome"
	
conda_env:
	@if [ -z "$$CONDA_DEFAULT_ENV" ]; then \
	echo "Please load bioinfo environment"; \
	exit 1; \
	fi

variables:
	@echo "ACC = $(ACC)"
	@echo "SPECIES = $(SPECIES)"
	@echo "N = $(N)"
	@echo "L = $(L)"
	@echo "R1 = $(R1)"
	@echo "R2 = $(R2)"
	@echo "SRA = $(SRA)"
	@echo "output_dir = $(output_dir)"
	@echo "read1 = $(read1)"
	@echo "read2 = $(read2)"
	@echo "trimmed_read1 = $(trimmed_read1)"
	@echo "trimmed_read2 = $(trimmed_read2)"

all: genome download trim index align call

clean: conda_env
	rm -rf $(output_dir)
	rm -f ncbi_dataset.zip
	rm -rf ncbi_dataset
	rm -f fastp*
	rm -f md5sum.txt
	rm -f README.md
	rm -f $(SPECIES)*
	rm -rf vcf
	rm -rf results
	rm -rf idx
	rm -rf src 
	rm -rf snpeff.config

## DOWNLOAD ##

genome: conda_env
	echo "Downloading genome"; \
	datasets download genome accession $(ACC) --include genome,gff3; \
	unzip -o ncbi_dataset.zip; \
	ln -sf ncbi_dataset/data/$(ACC)/*.fna $(SPECIES).fna; \
	ln -sf ncbi_dataset/data/$(ACC)/*.gff $(SPECIES).gff; \
	echo "Generating summary stats for the genome"
	wc -c $(SPECIES).fna
	seqkit stats $(SPECIES).fna
	echo "Number of contigs:"
	seqkit seq -i $(SPECIES).fna | grep ">" | wc -l
	echo "Summary of contigs:"
	seqkit fx2tab -n -l $(SPECIES).fna | awk '{print $$1 "\t" $$NF}'
	samtools faidx $(SPECIES).fna

download: conda_env
	@echo "Downloading from $(SRA) to $(output_dir) directory"
	fastq-dump -X 1000 -F --outdir $(output_dir) --split-files $(SRA)
	@echo "Running fastqc on $(output_dir)/$(read1) and $(output_dir)/$(read2)"
	mkdir -p $(output_dir)/fastqc
	fastqc $(output_dir)/$(read1) $(output_dir)/$(read2) -o $(output_dir)/fastqc

## QC ##

trim: conda_env
	@echo "Trimming adapters and low quality bases"
	mkdir -p $(output_dir)/trimmed_fastq
	fastp --thread 4 \
	--in1 "$(output_dir)/$(read1)" \
	--in2 "$(output_dir)/$(read2)" \
	--out1 "$(output_dir)/trimmed_fastq/$(read1)_trimmed.fq.gz" \
	--out2 "$(output_dir)/trimmed_fastq/$(read2)_trimmed.fq.gz" \
	-g \
	-q 20 \
	-u 20
	@echo "Running fastqc on trimmed reads"
	mkdir -p $(output_dir)/trimmed_fastq/fastqc
	fastqc $(output_dir)/trimmed_fastq/$(read1)_trimmed.fq.gz \
	$(output_dir)/trimmed_fastq/$(read2)_trimmed.fq.gz \
	-o $(output_dir)/trimmed_fastq/fastqc

flagstat: conda_env align
	@echo "Running samtools flagstat"
	samtools flagstat $(output_dir)/aligned_bam/aligned_reads_sorted.bam > $(output_dir)/aligned_bam/aligned_reads_sorted_stats.txt

filter: conda_env genome download trim align flagstat
	@echo "Filtering reads"
	mkdir -p $(output_dir)/filtered_bam
	samtools view -b -q 10 $(output_dir)/aligned_bam/aligned_reads_sorted.bam > $(output_dir)/filtered_bam/filtered_reads.bam
	@echo "Running samtools flagstat on filtered reads"
	samtools flagstat $(output_dir)/filtered_bam/filtered_reads.bam > $(output_dir)/filtered_bam/filtered_reads_stats.txt

## ALIGN ##

index: conda_env
	@echo "Indexing the genome"
	bwa index $(SPECIES).fna
	@echo "Generating the genome dictionary"
	picard CreateSequenceDictionary R=$(SPECIES).fna O=$(SPECIES).dict

align: conda_env genome download index trim
	@echo "Aligning reads to the genome"
	mkdir -p $(output_dir)/aligned_bam
	bwa mem -t 4 $(SPECIES).fna $(output_dir)/trimmed_fastq/$(trimmed_read1) $(output_dir)/trimmed_fastq/$(trimmed_read2) | samtools view -b - > $(output_dir)/aligned_bam/aligned_reads.bam
	@echo "Sorting the aligned reads"
	samtools sort -o $(output_dir)/aligned_bam/aligned_reads_sorted.bam $(output_dir)/aligned_bam/aligned_reads.bam
	@echo "Indexing the sorted aligned reads"
	samtools index $(output_dir)/aligned_bam/aligned_reads_sorted.bam
	@echo "Running samtools flagstat"
	samtools flagstat $(output_dir)/aligned_bam/aligned_reads_sorted.bam > $(output_dir)/aligned_bam/aligned_reads_sorted_stats.txt

## SIMULATE ##

simulate: conda_env
	@echo "Simulating reads for the genome"
	mkdir -p reads
	wgsim -N $(N) -1 $(L) -2 $(L) -r 0 -R 0 -X 0 $(SPECIES).fna $(R1) $(R2)
	seqkit stats $(R1) $(R2)

simalign: conda_env simulate
	@echo "Aligning simulated reads to the genome"
	
	mkdir -p $(output_dir)/simulated_aligned_bam

	bwa mem -t 4 $(SPECIES).fna $(R1) $(R2) | samtools view -b - > $(output_dir)/simulated_aligned_bam/simulated_aligned_reads.bam

	@echo "Sorting the simulated aligned reads"
	samtools sort -o $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam $(output_dir)/simulated_aligned_bam/simulated_aligned_reads.bam

	@echo "Indexing the sorted simulated aligned reads"
	samtools index $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam

	@echo "Running samtools flagstat"
	samtools flagstat $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam > $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted_stats.txt

call:
	mkdir -p vcf
	bcftools mpileup -Ou -f $(SPECIES).fna $(output_dir)/aligned_bam/aligned_reads_sorted.bam | \
	bcftools call -mv -Ou -o vcf/$(SPECIES).vcf

${SNPEFF_DB}: ${GFF} ${REF}	
	mkdir -p ${IDX}/${LABEL}

	# Copy the files to the snpEff folder.
	cp -f ${REF} ${IDX}/${LABEL}/sequences.fa
	cp -f ${GFF} ${IDX}/${LABEL}/genes.gff

	# Make the configuration file.
	echo "${LABEL}.genome : ${LABEL}" >	snpeff.config

	# Build the database.
	snpEff build -dataDir ${IDX} -v ${LABEL}

build: ${SNPEFF_DB}
	@ls -lh ${SNPEFF_DB}

snpeff: ${SNPEFF_DB}
	mkdir -p results
	snpEff ann -csvStats ${CSV} -s ${HTML} -dataDir ${IDX} -v ${LABEL} ${VCF} | bcftools view -O z -o ${EFF}
	bcftools index ${EFF}
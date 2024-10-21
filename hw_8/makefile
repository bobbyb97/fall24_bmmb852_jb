# Variables
ACC = GCA_000287275.1
SPECIES = Carsonella_ruddii
N = 1000
L = 100
R1 = reads/wgsim_read1.fq
R2 = reads/wgsim_read2.fq
SRA = DRR591238
output_dir = reads
read1 = $(SRA)_1.fastq
read2 = $(SRA)_2.fastq
trimmed_read1 = $(read1)_trimmed.fq.gz
trimmed_read2 = $(read2)_trimmed.fq.gz

# Targets
.PHONY: usage conda_env genome simulate download trim clean all

usage:
	@echo "usage: make [target]"
	@echo "Targets:"
	@echo "  usage    - Print the help"
	@echo "  conda_env- Check if bioinfo environment is loaded"
	@echo "  genome   - Download the genome"
	@echo "  simulate - Simulate reads for the genome"
	@echo "  download - Download reads from SRA"
	@echo "  trim     - Trim reads"
	@echo "  clean    - Remove all files"
	@echo "  all      - Run all targets"
	

conda_env:
	@if [ -z "$$CONDA_DEFAULT_ENV" ]; then \
	echo "Please load bioinfo environment"; \
	exit 1; \
	fi

genome: conda_env
	echo "Downloading genome"; \
	datasets download genome accession $(ACC); \
	unzip -o ncbi_dataset.zip; \
	ln -sf ncbi_dataset/data/$(ACC)/*.fna $(SPECIES).fna; \
	echo "Generating summary stats for the genome"
	wc -c $(SPECIES).fna
	seqkit stats $(SPECIES).fna
	echo "Number of contigs:"
	seqkit seq -i $(SPECIES).fna | grep ">" | wc -l
	echo "Summary of contigs:"
	seqkit fx2tab -n -l $(SPECIES).fna | awk '{print $$1 "\t" $$NF}'
	samtools faidx $(SPECIES).fna

simulate: conda_env
	@echo "Simulating reads for the genome"
	mkdir -p reads
	wgsim -N $(N) -1 $(L) -2 $(L) -r 0 -R 0 -X 0 $(SPECIES).fna $(R1) $(R2)
	seqkit stats $(R1) $(R2)

download: conda_env
	@echo "Downloading from $(SRA) to $(output_dir) directory"
	fastq-dump -X 1000 -F --outdir $(output_dir) --split-files $(SRA)
	@echo "Running fastqc on $(output_dir)/$(read1) and $(output_dir)/$(read2)"
	mkdir -p $(output_dir)/fastqc
	fastqc $(output_dir)/$(read1) $(output_dir)/$(read2) -o $(output_dir)/fastqc

trim: conda_env
	@echo "Trimming adapters and low quality bases"
	mkdir -p $(output_dir)/trimmed_fastq
	fastp --thread 4 \
	--in1 "$(output_dir)/$(read1)" \
	--in2 "$(output_dir)/$(read2)" \
	--out1 "$(output_dir)/trimmed_fastq/$(read1)_trimmed.fq.gz" \
	--out2 "$(output_dir)/trimmed_fastq/$(read2)_trimmed.fq.gz" \
	-g \
	-q 35 \
	-u 20
	@echo "Running fastqc on trimmed reads"
	mkdir -p $(output_dir)/trimmed_fastq/fastqc
	fastqc $(output_dir)/trimmed_fastq/$(read1)_trimmed.fq.gz \
	$(output_dir)/trimmed_fastq/$(read2)_trimmed.fq.gz \
	-o $(output_dir)/trimmed_fastq/fastqc


index: conda_env
	@echo "Indexing the genome"
	bwa index $(SPECIES).fna
	@echo "Generating the genome dictionary"
	picard CreateSequenceDictionary R=$(SPECIES).fna O=$(SPECIES).dict

align: conda_env
	@echo "Aligning reads to the genome"
	mkdir -p $(output_dir)/aligned_bam
	bwa mem -t 4 $(SPECIES).fna $(output_dir)/trimmed_fastq/$(trimmed_read1) $(output_dir)/trimmed_fastq/$(trimmed_read2) | samtools view -b - > $(output_dir)/aligned_bam/aligned_reads.bam
	@echo "Sorting the aligned reads"
	samtools sort -o $(output_dir)/aligned_bam/aligned_reads_sorted.bam $(output_dir)/aligned_bam/aligned_reads.bam
	@echo "Indexing the sorted aligned reads"
	samtools index $(output_dir)/aligned_bam/aligned_reads_sorted.bam
	@echo "Running samtools flagstat"
	samtools flagstat $(output_dir)/aligned_bam/aligned_reads_sorted.bam > $(output_dir)/aligned_bam/aligned_reads_sorted_stats.txt

simalign: conda_env
	@echo "Aligning simulated reads to the genome"
	mkdir -p $(output_dir)/simulated_aligned_bam
	bwa mem -t 4 $(SPECIES).fna $(R1) $(R2) | samtools view -b - > $(output_dir)/simulated_aligned_bam/simulated_aligned_reads.bam
	@echo "Sorting the simulated aligned reads"
	samtools sort -o $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam $(output_dir)/simulated_aligned_bam/simulated_aligned_reads.bam
	@echo "Indexing the sorted simulated aligned reads"
	samtools index $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam
	@echo "Running samtools flagstat"
	samtools flagstat $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted.bam > $(output_dir)/simulated_aligned_bam/simulated_aligned_reads_sorted_stats.txt

all: genome simulate download trim index align simalign

clean: conda_env
	rm -rf $(output_dir)
	rm -f ncbi_dataset.zip
	rm -rf ncbi_dataset
	rm -f fastp*
	rm -f md5sum.txt
	rm -f README.md
	rm -f $(SPECIES)*
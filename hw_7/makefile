# Variables
ACC = GCA_041682495.1
SPECIES = B_fervidus
N = 1000
L = 100
R1 = reads/wgsim_read1.fq
R2 = reads/wgsim_read2.fq
SRA = SRR28227798
output_dir = reads
read1 = $(SRA)_1.fastq
read2 = $(SRA)_2.fastq

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

all: genome simulate download trim

clean: conda_env
	rm -rf $(output_dir)
	rm -f ncbi_dataset.zip my_genome.gff3
	rm -f $(SPECIES).fna
	rm -rf ncbi_dataset
	rm -f fastp*
	rm -f md5sum.txt
	rm -f README.md
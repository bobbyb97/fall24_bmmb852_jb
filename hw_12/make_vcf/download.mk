-include config.mk

${REF}:
	mkdir -p ${DOWNLOAD_DIR}
	@echo "Downloading genome"
	micromamba run -n bioinfo datasets download genome accession ${ACC} --include genome,gff3 --filename ${DOWNLOAD_DIR}/${SPECIES}.zip
	unzip -o ${DOWNLOAD_DIR}/${SPECIES}.zip -d ${DOWNLOAD_DIR}/${SPECIES}_dataset
	cp ${DOWNLOAD_DIR}/${SPECIES}_dataset/ncbi_dataset/data/${ACC}/*.fna ${DOWNLOAD_DIR}/${SPECIES}.fna
	cp ${DOWNLOAD_DIR}/${SPECIES}_dataset/ncbi_dataset/data/${ACC}/*.gff ${DOWNLOAD_DIR}/${SPECIES}.gff

${DOWNLOAD_DIR}/${SPECIES}.stats: ${REF}
	micromamba run -n bioinfo seqkit stats ${DOWNLOAD_DIR}/${SPECIES}.fna > ${DOWNLOAD_DIR}/${SPECIES}.stats
	@echo "Number of contigs:" >> ${DOWNLOAD_DIR}/${SPECIES}.stats
	micromamba run -n bioinfo seqkit seq -i ${DOWNLOAD_DIR}/${SPECIES}.fna | grep ">" | wc -l >> ${DOWNLOAD_DIR}/${SPECIES}.stats
	@echo "Summary of contigs:" >> ${DOWNLOAD_DIR}/${SPECIES}.stats
	micromamba run -n bioinfo seqkit fx2tab -n -l ${DOWNLOAD_DIR}/${SPECIES}.fna | awk '{print $$1 "\t" $$NF}' >> ${DOWNLOAD_DIR}/${SPECIES}.stats

${SRA_DIR}:
	@echo "Downloading from ${SRA} to ${SRA_DIR}"
	micromamba run -n bioinfo fastq-dump --gzip -F --outdir ${SRA_DIR} --split-files ${SRA}
	@echo "Running fastqc on ${SRA_DIR}/${read1} and ${SRA_DIR}/${read2}"
	mkdir -p ${SRA_DIR}/fastqc
	micromamba run -n bioinfo fastqc ${SRA_DIR}/${read1} ${SRA_DIR}/${read2} -o ${SRA_DIR}/fastqc

all: ${REF} ${DOWNLOAD_DIR}/${SPECIES}.stats ${SRA_DIR}
	@echo "All done with downloads"
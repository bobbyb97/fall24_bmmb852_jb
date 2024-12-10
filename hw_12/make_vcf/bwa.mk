-include config.mk 
## ALIGN ##
${BAM}: ${TRIM_DIR}
	@echo "Indexing the genome"
	micromamba run -n bioinfo bwa index ${REF}

	@echo "Generating the genome dictionary"
	micromamba run -n bioinfo picard CreateSequenceDictionary R=${REF} O=${DOWNLOAD_DIR}/${SPECIES}.dict

	@echo "Aligning reads to the genome"
	mkdir -p ${BAM_DIR}
	micromamba run -n bioinfo bwa mem -t 4 ${REF} ${TRIM_DIR}/${trimmed_read1} ${TRIM_DIR}/${trimmed_read2} | micromamba run -n bioinfo samtools view -b > ${BAM}


${FLAGSTAT}: ${BAM}
	@echo "Running samtools flagstat"
	micromamba run -n bioinfo samtools flagstat ${BAM} > ${FLAGSTAT}

all: ${BAM} ${FLAGSTAT}
	@echo "All done with alignment"
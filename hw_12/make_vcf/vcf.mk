-include config.mk

${VCF}:
	mkdir -p ${VCF_DIR}
	micromamba run -n bioinfo bcftools mpileup -Ou -f ${REF} ${BAM} | \
	micromamba run -n bioinfo bcftools call -mv -Ou -o ${VCF_DIR}/${SPECIES}.vcf

${SNPEFF_DB}: ${GFF} ${REF}	
	mkdir -p ${IDX}/${LABEL}

	# Copy the files to the snpEff folder.
	cp -f ${REF} ${IDX}/${LABEL}/sequences.fa
	cp -f ${GFF} ${IDX}/${LABEL}/genes.gff

	# Make the configuration file.
	echo "${LABEL}.genome : ${LABEL}" >	${VCF_DIR}/snpeff.config

	# Build the database.
	micromamba run -n bioinfo snpEff build -dataDir ${IDX} -v ${LABEL}

build: ${SNPEFF_DB}
	@ls -lh ${SNPEFF_DB}

${EFF}: ${SNPEFF_DB} ${VCF}
	mkdir -p results
	micromamba run -n bioinfo snpEff ann -csvStats ${CSV} -s ${HTML} -dataDir ${IDX} -v ${LABEL} ${VCF} | micromamba run -n bioinfo bcftools view -O z -o ${EFF}
	micromamba run -n bioinfo bcftools index ${EFF}

all: ${EFF}
	@echo "All done with variant calling"
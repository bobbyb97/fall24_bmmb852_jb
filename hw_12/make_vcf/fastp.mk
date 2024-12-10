-include config.mk

${TRIM_DIR}: ${SRA_DIR}
	mkdir -p ${TRIM_DIR}

	@echo "Trimming adapters and low quality bases"
	micromamba run -n bioinfo fastp --thread 4 \
	--in1 "${SRA_DIR}/${read1}" \
	--in2 "${SRA_DIR}/${read2}" \
	--out1 "${TRIM_DIR}/${trimmed_read1}" \
	--out2 "${TRIM_DIR}/${trimmed_read2}" \
	--trim_poly_g \
	-q 20 \
	-u 20 \
	-h ${TRIM_DIR}/${SRA}_trimmed.html \
	-j ${TRIM_DIR}/${SRA}_trimmed.json

	@echo "Running fastqc on trimmed reads"
	mkdir -p ${TRIM_DIR}/fastqc
	micromamba run -n bioinfo fastqc ${TRIM_DIR}/${trimmed_read1} \
	${TRIM_DIR}/${trimmed_read2} \
	-o ${TRIM_DIR}/fastqc

all: ${TRIM_DIR}
	@echo "All done with trimming"
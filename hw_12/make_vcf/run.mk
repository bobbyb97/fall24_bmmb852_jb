-include config.mk

all:
	make -f ${SCRIPT_DIR}/download.mk all
	make -f ${SCRIPT_DIR}/fastp.mk all
	make -f ${SCRIPT_DIR}/bwa.mk all
	make -f ${SCRIPT_DIR}/index.mk all
	make -f ${SCRIPT_DIR}/vcf.mk all

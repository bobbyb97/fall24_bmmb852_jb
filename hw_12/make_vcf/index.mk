-include config.mk

${SORTED_BAM}: ${BAM}
	@echo "Sorting the aligned reads"
	samtools sort -o ${SORTED_BAM} ${BAM}
	@echo "Indexing the sorted aligned reads"
	samtools index ${SORTED_BAM}

all: ${SORTED_BAM}
	@echo "All done with sorting and indexing"
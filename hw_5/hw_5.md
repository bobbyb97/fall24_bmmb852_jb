## Generating Simulated FASTQ reads
-----
##### 1. Select a genome, then download the corresponding FASTA file.

I selected the genome of a virus, Vairimorpha bombi.
```
ACC="GCA_036780745.1"
SPECIES="V_bombi"

datasets download genome accession  
unzip -o ncbi_dataset.zip
ln -sf ncbi_dataset/data/${ACC}/*.fna ${SPECIES}.fna
```

The size of the file
    
```
wc -c ${SPECIES}.fna
```
* File size: 4803679 or 4.8 MB

The total size of the genome and number of contigs in the genome:
```
seqkit stats ${SPECIES}.fna

echo "Number of contigs:"
seqkit seq -i ${SPECIES}.fna | grep ">" | wc -l
```

* Genome length: 4,738,936

* Number of contigs: 57

The ID and length of each contig in the genome:
```
seqkit fx2tab -n -l ${SPECIES}.fna | awk '{print $1 "\t" $NF}' | head
```
Output:
```
JAWUGJ010000001.1       581600
JAWUGJ010000002.1       382231
JAWUGJ010000003.1       369002
JAWUGJ010000004.1       320259
JAWUGJ010000005.1       274745
JAWUGJ010000006.1       245381
JAWUGJ010000007.1       215746
JAWUGJ010000008.1       180348
JAWUGJ010000009.1       176623
JAWUGJ010000010.1       171846
```

##### 2. Generate a simulated FASTQ output for a sequencing instrument of your choice.  Set the parameters so that your target coverage is 10x.

Generating simulated FASTQ using wgsim
```
wgsim -N ${N} -1 ${L} -2 ${L} -r 0 -R 0 -X 0 ${SPECIES}.fna ${R1} ${R2}
```

How many reads have you generated?

* 476644

What is the average read length?

* 100

How big are the FASTQ files?

* 58 MB each

Compress the files and report how much space that saves.

* Each file is 10 MB now, so in total I saved 96 MB
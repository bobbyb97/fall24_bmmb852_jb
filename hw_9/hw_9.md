#### HW 9: Filtering a BAM file
----
1.  How many reads did not align with the reference genome?
```
samtools view -c -f 4 reads/aligned_bam/aligned_reads.bam
```
Output:
```
1304
```
2. How many primary, secondary, and supplementary alignments are in the BAM file?
```
# Primary 
samtools view -c -F 256 reads/aligned_bam/aligned_reads.bam 
# Secondary
samtools view -c -f 256 reads/aligned_bam/aligned_reads.bam
# Supplementary
samtools view -c -f 2048 reads/aligned_bam/aligned_reads.bam
```

Output:
```
Primary: 1979

Secondary: 0

Supplementary: 1
```


3. How many properly-paired alignments on the reverse strand are formed by reads contained in the first pair ?
```
samtools view -c -f 2 reads/aligned_bam/aligned_reads.bam
```
Output:
```
625
```

Make a new BAM file that contains only the properly paired primary alignments with a mapping quality of over 10
```
samtools view -b -q 10 $(output_dir)/aligned_bam/aligned_reads_sorted.bam > $(output_dir)/filtered_bam/filtered_reads.bam
```

Original flagstats:
```
1979 + 0 in total (QC-passed reads + QC-failed reads)
1978 + 0 primary
0 + 0 secondary
1 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
675 + 0 mapped (34.11% : N/A)
674 + 0 primary mapped (34.07% : N/A)
1978 + 0 paired in sequencing
989 + 0 read1
989 + 0 read2
624 + 0 properly paired (31.55% : N/A)
626 + 0 with itself and mate mapped
48 + 0 singletons (2.43% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```

Filtered flagstats:
```
629 + 0 in total (QC-passed reads + QC-failed reads)
628 + 0 primary
0 + 0 secondary
1 + 0 supplementary
0 + 0 duplicates
0 + 0 primary duplicates
629 + 0 mapped (100.00% : N/A)
628 + 0 primary mapped (100.00% : N/A)
628 + 0 paired in sequencing
316 + 0 read1
312 + 0 read2
578 + 0 properly paired (92.04% : N/A)
580 + 0 with itself and mate mapped
48 + 0 singletons (7.64% : N/A)
0 + 0 with mate mapped to a different chr
0 + 0 with mate mapped to a different chr (mapQ>=5)
```
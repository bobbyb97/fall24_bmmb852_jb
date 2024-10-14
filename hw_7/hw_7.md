## HW 7: Makefiles
---
The makefile contains 8 targets, all aimed towards downloading a genome from NCBI, simulating reads, downloading FASTQ reads, and quality control.

### Targets:
#### `usage`
Prints the help message, listing all available targets and their descriptions.

#### `conda_env`
Checks to see if the bioinfo env is loaded. Each following target depends on this.

#### `genome` 
Downloads the genome specified by the $ACC variable and saves it under the specified $SPECIES name.

#### `simulate` 
Simulates reads ($R1, $R2) based on number ($N), length ($L), and saves them to /reads/.

#### `download` 
Downloads FASTQ files from $SRA and saves to $output_dir.

#### `trim`
Trims FASTQ files using fastp.

#### `clean` 
Removes all files generated using any of the commands in the makefile.

#### `all`
Runs all targets except clean.
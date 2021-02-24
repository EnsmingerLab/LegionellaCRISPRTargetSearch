# Legionella CRISPR Target Search

# Introduction

Legionella CRISPR Target Search is a shell script that was designed to search for targets of Legionella pneumophila's CRISPR-Cas systems in genomic data.

It runs on any Unix operating system with the software dependancies mentioned below.

# Dependancies

The following packages need to be installed for the script to run:

| Dependancy | Link for download
--- | --- |
| blastn | Comes with ncbi-blast+ package. Download from https://ftp.ncbi.nlm.nih.gov/blast/executables/blast+/LATEST/
| bedtools | Download from https://github.com/arq5x/bedtools2
| samtools | Download from https://github.com/samtools/samtools
| GNU parallel | Download from https://www.gnu.org/software/parallel/


Alternatively, you can use the homebrew package manager for MacOS to install the dependancies.

```
brew install blast bedtools samtools parallel
```

# Basic command

```
bash Legionella_CRISPR_Target_Search.sh <SPACER LIST HERE> <ALIGNMENT LENGTH FILTER HERE> <MISMATCH FILTER HERE>
```

### Parameters

Spacer list:	Ensure the spacer list is a FASTA formatted file with a .txt extension. Ensure it's in the same directory as the script file and the input genomes.

Alignment length filter:	A numerical value that will be used to determine an alignment length cutoff when filtering BLAST results.

Mismatch filter:	A numerical value that will be used to determine a mismatch cutoff when filtering BLAST results.

# Additional information

### Link to bioRxiv submission that uses this script

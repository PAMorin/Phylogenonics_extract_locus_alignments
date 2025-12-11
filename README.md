# Phylogenomics extract locus alignments

This README describes a pipeline for identifying a set of gene (or other locus) sequences in a reference genome to map their locations, then using the locations to extract the sequences from a set of bam files and create a set of locus alignments for phylogenomic analysis.

## Description



## Getting Started

### Dependencies

Before running the pipeline, make sure you have installed the following programs: 

bwa-mem2/2.2.1
samtools/1.19


### Installing

If you want to run this pipeline, you can simply clone the github repository using the code below:
```
git clone git@github.com:pmorin?????
```

## Executing program

### Step 1: map locus sequences to a reference genome and generate a bed file for their coordinates.
UCE_map2ref_extract_bed.sh

**Required inputs:**
Fasta file of locus sequences (example = uce-5k-probes.fasta)

**Outputs**
Six-column bed file of locus coordinates (including optional sequence padding on each end of the locus coordinates)
### output bed file has 6 columns:
1. chromosome ID
2. start position
3. end position
4. locus name
5. mapping quality
6. strand (forward = '+', reverse = '-')

**IMPORTANT**
In the script, change the paths, pad, and output file name:
REFERENCE_GENOME="path_to_your_reference_genome";
QUERY_SEQUENCES="path_to_your_query_sequences";
OUTPUT_DIR="path_to_your_output_directory";
PAD=(int);
OUTPUT_BED="your_bed_file_name";

### Step 2

**Required inputs:**

**Outputs**

**IMPORTANT**

### Step 3

**Required inputs:**

**Outputs**

**IMPORTANT**

### Step 4

**Required inputs:**

**Outputs**

**IMPORTANT**



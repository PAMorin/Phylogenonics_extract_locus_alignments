# Phylogenomics extract locus alignments

This README describes a pipeline for identifying a set of gene (or other locus) sequences in a reference genome to map their locations, then using the locations to extract the sequences from a set of bam files and create a set of locus alignments for phylogenomic analysis.

## Description



## Getting Started

### Dependencies

Before running the pipeline, make sure you have installed the following programs: 

bwa-mem2/2.2.1  
samtools/1.19  
bedtools/2.31.1  
htslib/1.19 (part of SAMtools and BCFtools, but may need to be installed separately to use bgzip, htsfile, and tabix utilities)


### Installing

To copy this pipeline from the github repository https://github.com/PAMorin/Phylogenomics_extract_locus_alignments:

1. Click the green "< > Code" button above the file list.
2. Select Download ZIP from the dropdown menu. 


## Executing program

### Step 1: Map locus sequences to a reference genome and generate a bed file for their coordinates.
```
sbatch map2ref_extract_bed.sh
```
**Required inputs:**  
Fasta file of locus sequences (example = uce-5k-probes.fasta)

**Outputs**  
Six-column bed file of locus coordinates (including optional sequence padding on each end of the locus coordinates)
### Output bed file has 6 columns:
1. chromosome ID
2. start position
3. end position
4. locus name
5. mapping quality
6. strand (forward = '+', reverse = '-')

**IMPORTANT**  
In the script, change the paths, pad, and output file name:
REFERENCE_GENOME="path_to_your_reference_genome"  
QUERY_SEQUENCES="path_to_your_query_sequences"  
OUTPUT_DIR="path_to_your_output_directory"  
PAD=(int)  
OUTPUT_BED="your_bed_file_name"  


### Step 2: generate consensus genomes from mpileup (bam) files for each sample
```
sbatch consensus_genome_woIUPAC.sh
```
**Required inputs:**  
Bam file for each species whole genome sequence data aligned to the reference genome  
Indexed reference genome  

**Outputs**  
Consensus genome sequences  
(Each nucleotide call based on the most common base at the site. In the case of ties a random base is chosen among the bases with the same maximum counts. N's or filtered based are ignored. If multiple individuals are used the four bases are counted across individuals. Indels are ignored, so that the consensus sequence is the same length as the reference sequence.)

**IMPORTANT**  
The default organization and naming structure required for the script is for each bam file name to start with <species_sample>, and to be in it's own subdirectory that is also named <species_sample> (e.g., the path to the bam file for <Felis_catus> would be /maindir/Felis_catus/Felis_catus_dedup.bam).

### Step 3: Extract the locus sequences from each consensus sequence
```
sbatch extract_loci_array.sh
```
**Required inputs:**  
Consensus sequence from Step 2  
Bed file from Step 1

**Outputs**  
Fasta file for each sample, containing all of the extracted loci.


### Step 4: Extract each locus from each sample file and combine them into a locus alignment of all samples in individual locus files. 

**Required inputs:**  
Text file including just the "description" column from the bed file (e.g., UCE_loci.txt)  
Directory of fasta files containing extracted loci for each sample from Step 3

**Outputs**  
Directory of locus alignments for each locus, containing the locus sequences for all samples.

### Step 5: Filter files to remove locus alignments with more than a specified portion of N's in the alignment.
```
sbatch Filter_XpctN_alignments.sh
```
**Required inputs:**  
Directory of locus alignments for each locus from Step 4.

**Outputs**  
Two directories containing the locus alignment with greater-than and less-than the specified percent N cutoff. 



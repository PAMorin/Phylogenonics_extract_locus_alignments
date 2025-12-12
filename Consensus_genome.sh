#!/bin/bash

#SBATCH --job-name=consensus_genome   ## job name
#SBATCH -e consensus_genome%j.e.txt    ## error message name
#SBATCH -o consensus_genome.log.%j.out  ## log file output name
#SBATCH -c 10    ## <number of cores to ask for> (cpu cores per task?)
#SBATCH -p medmem
#SBATCH --mem=80G
#SBATCH --array=1-2
#SBATCH -t 48:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 

##########################################################################################
# Uses bcftools to extract consensus genome from bam files in parallel (Slurm array). All bam files are in separate  subdirectories within a main directory, named as "species_sample". Bam files are named as "species_sample_dedup.bam", so that "species_sample" is the same for the directory and the first part of the file name. 

module load bio/samtools/1.19
module load bio/bcftools/1.11
module load tools/rclone
module load tools/pigz/2.4

set -eux

#### --- User-defined variables ---
# the BAMLIST is a list of directories, and assumes each bam directory is in a single directory, and each bam file is in its own subdirectory with the same name scheme as the files.  
BAMLIST=(
B.b.acu_DRR014695
B.b.bor_SRR26062094
)

# reference genome used for WGS read alignments
REFDIR=/home/pmorin/Ref_genomes/Egla/GCA_028564815.1_mEubGla1_pri
REF=GCA_028564815.1_mEubGla1.hap2_genomic.fna
NGS=/home/pmorin/projects/Miscellaneous/mysticete_phylo/fastq_files
CONS_OUT=/home/pmorin/projects/Miscellaneous/TEST_phylogenomics_extract_loc_align_repo/consensus_genomes
baseQ=30
mapQ=25
bam=dedup.bam # supply appropriate wildcard to find the bam file within BAMDIR (e.g., if your bam files in each directory end with "dedup.bam", then files with *dedup.bam ending will be used to generate the consensus)
THREADS=10

#### --- end User-defined variables ---

# make output directories if they don't already exist
mkdir -p ${CONS_OUT}

#####################

# generate windows list for each scaffold
BAMDIR=${BAMLIST[$SLURM_ARRAY_TASK_ID-1]} # if using an alias list within the script
echo ${BAMDIR}

TEMP=${CONS_OUT}/${BAMDIR}/TEMP
mkdir -p $TEMP

# extract the species and sample ID's from the bam directory, assuming name format = species_sample.
SP=`echo $BAMDIR | cut -f1 -d "_"`
sample=`echo $BAMDIR | cut -f2 -d "_"`
BAMFILE=$(ls ${NGS}/${BAMDIR}/*${bam}) 
echo ${BAMFILE}
BAMNAME=$(basename "$BAMFILE")
echo ${BAMNAME}

# Define output file name prefix
OUT_PREFIX=${SP}_${sample}_consensus_genome

#########
# call variants
# based on http://samtools.github.io/bcftools/howtos/consensus-sequence.html and https://www.biostars.org/p/353780/#485121

# 1. Combine pileup, calling, and normalization into one stream, writing to a file
bcftools mpileup -Ou -f ${REFDIR}/${REF} ${BAMFILE} | \
bcftools call --threads ${THREADS} -Ou -mv | \
bcftools norm --threads ${THREADS} -f ${REFDIR}/${REF} -Oz -o ${TEMP}/${BAMNAME}_output.vcf.gz

# 2. Index the resulting compressed VCF (Required for consensus)
bcftools index --threads ${THREADS} ${TEMP}/${BAMNAME}_output.vcf.gz

# 3. Generate the consensus sequence
bcftools consensus --iupac-codes -f ${REFDIR}/${REF} ${TEMP}/${BAMNAME}_output.vcf.gz > ${CONS_OUT}/${OUT_PREFIX}.fa

# delete vcf file to save space
rm -r ${TEMP}

# -I, --iupac-codes
# output variants in the form of IUPAC ambiguity codes determined from FORMAT/GT fields. By default all samples are used and can be subset with -s, --samples and -S, --samples-file. Use -s - to ignore samples and use only the REF and ALT columns. NOTE: prior to version 1.17 the IUPAC codes were determined solely from REF,ALT columns and sample genotypes were not considered.
# Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion
# bcftools mpileup -Ou -f ref.fa aln.bam collects metrics about each covered reference base. This information are used to find suspicious sites
# bcftools call -Ou -mv calls variants on suspicious sites
# bcftools norm -f ref.fa normalize variants
# bcftools CONS_OUT -f ref.fa -o CONS_OUT.fa creates a CONS_OUT.fa based on the variants found
##########

# compress fasta files
pigz -p ${THREADS} ${CONS_OUT}/${OUT_PREFIX}.fa

# index fasta files
samtools faidx ${CONS_OUT}/${OUT_PREFIX}.fa.gz



#!/bin/bash

#SBATCH --job-name=CONS_OUT_genome   ## job name
#SBATCH -e CONS_OUT_genome%j.e.txt    ## error message name
#SBATCH -o CONS_OUT_genome.log.%j.out  ## log file output name
#SBATCH -c 10    ## <number of cores to ask for> (cpu cores per task?)
#SBATCH -p medmem
#SBATCH --mem=80G
#SBATCH --array=1-2
#SBATCH -t 48:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 

##########################################################################################
# Uses Angsd to extract CONS_OUT genome

module load bio/samtools/1.19
module load bio/bcftools/1.11
module load tools/rclone
module load tools/pigz/2.4

set -eux

# --- User-defined variables ---
# the BAMLIST assumes each bam directory is in a single directory, and each bam file is in its own subdirectory with the same name scheme as the files 
BAMLIST=(
B.b.acu_DRR014695
B.b.bor_SRR26062094
)

REFDIR=/home/pmorin/Ref_genomes/Egla/GCA_028564815.1_mEubGla1_pri
REF=GCA_028564815.1_mEubGla1.hap2_genomic.fna
THREADS=10

NGS=/home/pmorin/projects/Miscellaneous/mysticete_phylo/fastq_files
CONS_OUT=/home/pmorin/projects/Miscellaneous/TEST_phylogenomics_extract_loc_align_repo/consensus_genomes
TEMP=${CONS_OUT}/TEMP
baseQ=30
mapQ=25
bam=dedup.bam # supply appropriate wildcard to find the bam file within BAMDIR (e.g., if your bam files in each directory end with "dedup.bam", then files with *dedup.bam ending will be used to generate the consensus)

# --- end User-defined variables ---

# make output directories if they don't already exist
mkdir -p ${CONS_OUT}
mkdir -p ${TEMP}

#####################

# Get working bam file based on array number
NUM=$(printf %02d ${SLURM_ARRAY_TASK_ID})

# generate windows list for each scaffold
BAMDIR=${BAMLIST[$SLURM_ARRAY_TASK_ID-1]} # if using an alias list within the script
echo ${BAMDIR}

# extract the species and sample ID's from the bam directory, assuming name format = Gspe_sample, for Species_sample.
SP=`echo $BAMDIR | cut -f1 -d "_"`
sample=`echo $BAMDIR | cut -f2 -d "_"`
BAMFILE=$(ls ${NGS}/${BAMDIR}/*${bam}) 
echo ${BAMFILE}
BAMNAME=`echo $BAMFILE | cut -f7 -d "/"` # alter to capture bam file name from path.

# Define output file name prefix
OUT_PREFIX=${SP}_${sample}_consensus_genome

#########
# call variants
# based on http://samtools.github.io/bcftools/howtos/consensus-sequence.html and https://www.biostars.org/p/353780/#485121

bcftools mpileup -Ou -f ${REFDIR}/${REF} ${BAMFILE} | bcftools call --threads ${THREADS} -Ou -mv | bcftools norm --threads ${THREADS} -f ${REF_GENOME} -Oz -o ${TEMP}/${BAMNAME}_output.vcf.gz
bcftools index --threads ${THREADS} ${TEMP}/${BAMNAME}_output.vcf.gz
bcftools consensus --iupac-codes -f ${REF_GENOME} ${TEMP}/${BAMNAME}_output.vcf.gz > ${CONS_OUT}/${OUT_PREFIX}.fa

# -I, --iupac-codes
# output variants in the form of IUPAC ambiguity codes determined from FORMAT/GT fields. By default all samples are used and can be subset with -s, --samples and -S, --samples-file. Use -s - to ignore samples and use only the REF and ALT columns. NOTE: prior to version 1.17 the IUPAC codes were determined solely from REF,ALT columns and sample genotypes were not considered.
# Use the -Ou option when piping between bcftools subcommands to speed up performance by removing unnecessary compression/decompression and VCF←→BCF conversion
# bcftools mpileup -Ou -f ref.fa aln.bam collects metrics about each covered reference base. This information are used to find suspicious sites
# bcftools call -Ou -mv calls variants on suspicious sites
# bcftools norm -f ref.fa normalize variants
# bcftools CONS_OUT -f ref.fa -o CONS_OUT.fa creates a CONS_OUT.fa based on the variants found
##########

# compress fasta files
pigz -p ${THREADS} *.fa

# index fasta files
samtools faidx *.fa.gz



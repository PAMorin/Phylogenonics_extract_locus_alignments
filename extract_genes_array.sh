#!/bin/bash
#SBATCH --job-name=extract_seqs     ## job name
#SBATCH -e extract_seqs_%j.e.txt    ## error message name
#SBATCH -o extract_seqs.log.%j.out  ## log file output name
#SBATCH -p medmem
#SBATCH -c 1
#SBATCH --mem=5G
#SBATCH --array=1-4 # %# limits the number of subjobs in an array that run concurrently, e.g., array 1-10, with 2 jobs at a time, would be --array=1-10%2.
#SBATCH -t 8:00:00   ## walltime in mins or mins:secs or hrs:mins:secs. 

####################################################################################

# Load necessary modules
module load bio/bedtools/2.31.1
module load bio/samtools/1.19
module load bio/htslib/1.19

# Define input and output directories and files
THREADS=5
FASTADIR=/home/pmorin/NGS_data_files/Delphinidae_phylo/Consensus_wo_iupac
UCEDIR=/home/pmorin/projects/Miscellaneous/Delphinidae_UCE
BEDFILE=$UCEDIR/UCE_Oorc_GCA_937001465.bed2_padded.txt
OUTPUTDIR=$UCEDIR/UCE_loci_out_stringent

# Consensus sequences
FASTAFILES=(
Ccom_SRR12437578_angsd_consensus_genome.fa.gz
Ccom_z480_angsd_consensus_genome.fa.gz
Ceut_z2317_angsd_consensus_genome.fa.gz
Chea_z7320_angsd_consensus_genome.fa.gz
)

mkdir -p $OUTPUTDIR

cd $FASTADIR

# Extract sequence based on array index
FASTA=${FASTAFILES[$SLURM_ARRAY_TASK_ID-1]}

FA=`echo $FASTA | cut -f1 -d "."` # fasta file name before ".fa.gz"
SP=`echo $FASTA | cut -f1 -d "_"`
SPID=`echo $FASTA | cut -f1,2 -d "_"`

echo "indexing $FASTA..."
# 
# # unzip fasta files # this shouldn't be necessary with bedtools 2.31.1 (fixed issue)
# # gunzip $FASTADIR/$FASTA
# 
# index fasta files
# samtools faidx $FASTADIR/$FA.fa.gz

echo "Extracting sequences from $FASTA..."

# Use bedtools to extract sequences from fasta file based on bed file coordinates
bedtools getfasta -fi $FASTADIR/$FA.fa.gz -bed $BEDFILE -fo $OUTPUTDIR/${FA%.fa}.bedout2.fa -name+
# -name+ adds information from additional columns to the fasta headers

# add species name to each fasta header
sed -E "s/^>/>${SPID}_/g" $OUTPUTDIR/${FA%.fa}.bedout2.fa > $OUTPUTDIR/${FA%.fa}.bedout.fa

# rm $OUTPUTDIR/${FA%.fa}.bedout2.fa # for some reason this causes some to fail, so leave off of script and remove manually


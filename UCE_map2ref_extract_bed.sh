#!/bin/bash

#SBATCH --job-name=UCE_map2ref   ## job name
#SBATCH -e UCE_map2ref_%j.e.txt    ## error message name
#SBATCH -o UCE_map2ref.log.%j.out  ## log file output name
#SBATCH -c 5
#SBATCH -p medmem
#SBATCH --mem=45G
#SBATCH -t 01:00:00   ## walltime in mins, mins:secs, hrs:mins:secs, days-hours 
#SBATCH --ntasks=1 
##########################################################################################
# map UCE sequences to a reference genome to find the sequence positions, then extract the positions and sequence IDs to a bed file.

# Load the desired aligners module
module load aligners/bwa-mem2/2.2.1
module load bio/samtools/1.19

set -eux # halt if error, missing value; print commands first prints the command to the terminal, using a + to indicate that the output is a command

# --- User-defined variables ---
# CHANGE THESE PATHS TO YOUR ACTUAL FILES AND DIRECTORIES
REFERENCE_GENOME="/home/pmorin/Ref_genomes/Oorc/NCBI_GCA_937001465.1/GCA_937001465.1_mOrcOrc1.1_genomic.fna"
QUERY_SEQUENCES="/home/pmorin/projects/Miscellaneous/TEST_phylogenomics_extract_loc_align_repo/uce-5k-probes.fasta"
OUTPUT_DIR="/home/pmorin/projects/Miscellaneous/TEST_phylogenomics_extract_loc_align_repo"
PAD=500  # Amount of base pairs to add to both sides of the alignment, e.g., for short loci like UCEs. Not needed for long loci (like full genes)
OUTPUT_BED="${OUTPUT_DIR}/UCE_Oorc_GCA_937001465.bed.txt" # Output file name for the BED data

# Ensure the output directory exists
mkdir -p ${OUTPUT_DIR}

# --- Indexing (only needed if not already indexed) ---
# Check if the bwa-mem2 index files exist. If not, run indexing.
if [ ! -f ${REFERENCE_GENOME}.bwt.2bit.64 ]; then
  echo "Indexing reference genome..."
  bwa-mem2 index ${REFERENCE_GENOME}
  echo "Indexing complete."
fi

# --- Alignment and Conversion ---
# 1. bwa-mem2: Maps sequences.
# 2. samtools view:
#    -F 4: Exclude unmapped.
#    -q 20: Exclude MAPQ < 20.
# 3. awk:
#    -v pad="$PAD": Imports the bash variable into awk.
#    - Parses CIGAR, determines Strand.
#    - Filters for length >= 100bp.
#    - Applies padding and prevents negative start coordinates.
# 4. sort: Sorts by Chromosome then Start position.

# output bed file has 6 columns:
# 1. chromosome ID
# 2. start position
# 3. end position
# 4. locus name
# 5. mapping quality
# 6. strand (forward = '+', reverse = '-')

bwa-mem2 mem -t $SLURM_CPUS_PER_TASK ${REFERENCE_GENOME} ${QUERY_SEQUENCES} | \
  samtools view -F 4 -q 20 - | \
  awk -v pad="$PAD" '
    {
      qname=$1;
      flag=$2;
      chr=$3;
      start=$4;
      mapq=$5;
      cigar=$6;
      end=start;

      # 1. Calculate End Position via CIGAR
      temp_cigar = cigar;
      while (match(temp_cigar, /([0-9]+)([MDN])/)) {
        len = substr(temp_cigar, RSTART, RLENGTH - 1);
        op = substr(temp_cigar, RSTART + RLENGTH - 1, 1);
        if (op == "M" || op == "D" || op == "N") {
          end += len;
        }
        temp_cigar = substr(temp_cigar, RSTART + RLENGTH);
      }

      # 2. Determine Strand
      if (int(flag / 16) % 2 == 1) {
        strand = "-";
      } else {
        strand = "+";
      }

      # 3. Filter, Pad, and Print
      # Only process if original alignment is >= 100bp
      if ((end - start) >= 100) {
        
        # Calculate standard BED coordinates (0-based start, 1-based end)
        bed_start = start - 1;
        bed_end = end - 1;

        # Apply padding using the passed variable
        pad_start = bed_start - pad;
        pad_end = bed_end + pad;

        # Safety Check: Ensure start is not negative
        if (pad_start < 0) {
          pad_start = 0;
        }

        # Output BED6: Chr, Padded_Start, Padded_End, Name, Score, Strand
        print chr, pad_start, pad_end, qname, mapq, strand
      }
    }
  ' OFS="\t" | \
  sort -k1,1 -k2,2n > ${OUTPUT_BED}
  
echo "Alignment and 6-column BED file creation complete. Output saved to: ${OUTPUT_BED}"
#!/bin/bash
#SBATCH --job-name=remove_pctN
#SBATCH -e remove_pctN_%j.e.txt    ## error message name
#SBATCH -o remove_pctN.log.%j.out  ## log file output name
#SBATCH --mail-user=phillip.morin@noaa.gov
#SBATCH --mail-type=ALL  ## (NONE, BEGIN, END, FAIL, ALL)
#SBATCH --ntasks-per-node=1
#  SBATCH --cpus-per-task=20
#SBATCH -c 1
#  SBATCH --mem-per-cpu=4G
#SBATCH --time=1:00:00
#SBATCH --output=fasta-processing.%j.out
#########################################################################################

# Filter the locus alignments to remove alignments with >X% N's, copying them to a "high-N" subdirectory, then copy those with ≤X% Ns to a "low_N" subdirectory. 

module load tools/rclone/1.59.2

#### --- User-defined variables ---

# Set path to directory containing the fasta files
FASTA_DIR=/home/pmorin/projects/Miscellaneous/TEST_phylogenomics_extract_loc_align_repo/aligned_loci_out

# Set names and path for output directories
highN_DIR=${FASTA_DIR}/loci_high_Ns
lowN_DIR=${FASTA_DIR}/loci_low_Ns

pctN=1 # percent N cutoff value

#### --- end user-defined variables ---

# make output directories
mkdir -p ${highN_DIR}
mkdir -p ${lowN_DIR}

# Change to the fasta directory
cd $FASTA_DIR

# Loop through all fasta files in the directory
for FILE in *.fasta
do
    # Get the total length of all sequences in the file
    LENGTH=$(grep -v '>' $FILE | tr -d '\n' | wc -c)
    
    # Get the number of N's in the file
    NUM_NS=$(grep -v '>' $FILE | tr -d '\n' | tr -cd 'Nn' | wc -c)
    
    # Calculate the percentage of N's in the file
    PERCENT_NS=$(echo "scale=4; ($NUM_NS/$LENGTH)*100" | bc)
    
    # If the percentage of N's is greater than pctN, rename the file
    if (( $(echo "$PERCENT_NS > $pctN" | bc -l) )); then
        mv $FILE ${highN_DIR}/${FILE}
    fi
done

mv $FASTA_DIR/*.fasta $lowN_DIR/

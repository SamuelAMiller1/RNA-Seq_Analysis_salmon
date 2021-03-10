#!/bin/bash

while [[ $# -gt 0 ]]; do
  key=$1
  case $key in
  
    -s|--seq-dir)
      seqdir="$2"
      shift 2
      ;;
    -a|--assembly)
      ASSEMBLY="$2"
      shift 2
      ;;
    -o|--out-dir)
      outdir="$2"
      shift 2
      ;;
    -t|--threads)
      CORES="$2"
      shift 2
      ;;

  esac
done


# Navigate to parent directory
cd $outdir

# Collect fastqs into array
FASTQS=($(find $seqdir -type f -name "*.fastq.gz"))

# Make results file

mkdir -p results/fastqc/

# Run FASTQC

fastqc \
 -t $CORES \
 -o results/fastqc \
 ${FASTQS[@]}

# Remove .zip files (if unnecessary to keep)

ZIP=($(find results/fastqc -type f -name "*fastqc.zip"))

for files in "${ZIP[@]}"; do
 rm -f "$files"
done


#########################
# INDEXING AND COUNTING #
#########################

#Create directories to store the index and quant matrix

mkdir -p genome_index
mkdir -p quants

##################################
## Generate Salmon Genome Index ##
##################################

# Index the Genome

salmon index -t $ASSEMBLY -i ./genome_index

##################################
## Generate Salmon Pseudocounts ##
##################################

# Quantify reads

# Create array for FASTQ files

FASTQS=($(find $seqdir -name "*\.fastq"))

#Subset FASTQ's in individual arrays for Read 1 and Read 2 

FASTQSONE=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep _R1_`)
FASTQSTWO=(`echo ${FASTQS[@]} | sed 's/ /\n/g' | grep _R2_`)


# Align the reads with STAR

for n in {0..${#FASTQSONE[@]}}; do

 R1=${FASTQSONE[${n}]}
 R2=${FASTQSTWO[${n}]}

 samp=`basename ${R1}`
 salmon quant -i ./genome_index -l A \
        -1 $R1 \
        -2 $R2 \
	-o quants/${samp} \
	--validateMappings \
	--seqBias
done




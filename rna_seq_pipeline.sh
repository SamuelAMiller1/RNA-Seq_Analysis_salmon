#!/bin/bash

# Activate enviornment

export PATH=/N/u/millesaa/Carbonate/miniconda3/envs/rnaseq/bin:$PATH

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

# Collect fastqs into array
FASTQS=($(find ./sequences -type f -name "*.fastq.gz"))

# Run FASTQ

fastq \
 -t 6 \
 ${FASTQS[@]}

# Make results file

mkdir -p results/fastqc/

# Move .html files

HTML=($(find ./sequences -type f -name "*fastqc.html"))

mv ${HTML[@]} results/fastqc/

# Remove .zip files

ZIP=($(find ./sequences -type f -name "*fastqc.zip"))

for files in "${ZIP[@]}"; do
 rm -f "$files"
done


_______________________




# Assign variables to annotation and assembly file locations

ASSEMBLY='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.transcripts.fa.gz'

ANNOTATION='ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_37/gencode.v37.annotation.gtf.gz'

# Navigate to parent directory

cd /N/slate/millesaa/dll1_corin_ht29_rnaseq/dll1_corin_rnaseq

#Create directories to store the index and quant matrix

mkdir -p genome
mkdir -p quants

##################################
## Generate Salmon Genome Index ##
##################################

# Download and unpack the genome assembly and annotation

curl $ASSEMBLY -o ./genome/gencode.v37.transcripts.fa.gz

# Index the Genome

salmon index -t ./genome/gencode.v37.transcripts.fa.gz -i ./genome/gencode_v37_index --gencode


##################################
## Generate Salmon Pseudocounts ##
##################################

# Quantify reads
# Next run, add flags --seqBias --useVBOpt

for fn in sequences/ILMN_906_Ohagan_{1..24};
 do
 samp=`basename ${fn}`
 salmon quant -i ./genome/gencode_v37_index -l A \
        -1 ${fn}/*\_R1_001.fastq.gz \
        -2 ${fn}/*\_R2_001.fastq.gz \
	-o quants/${samp} \
	--validateMappings \
	--seqBias
done


________________________________


#!/bin/bash

## Create a single .tsv file that contains gencode v37 transcript and gene ID's

# Copy transcript ID's to text file

zgrep -P -o 'ENST\d{11}' genome/gencode.v37.transcripts.fa.gz > enst.txt


# Copy gene ID's to text file

zgrep -P -o 'ENSG\d{11}' genome/gencode.v37.transcripts.fa.gz > ensg.txt

# Conjoin txt files into one tsv file

paste enst.txt ensg.txt > gene_map.tsv


________________________________________________



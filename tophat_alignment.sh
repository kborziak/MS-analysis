#!/bin/bash
#$ -S /bin/bash

# This script runs tophat on a batch of fastq.gz files
# paths to scripts and directories
TOPHAT = /home/borzik01/tophat-2.0.14.Linux_x86_64/tophat
BOWTIE_INDEX = home/borzik01/genomes/Homo_sapiens/NCBI_spiked/GRCh38/Sequence/Bowtie2Index/genome
GENOME = home/borzik01/genomes/Homo_sapiens/NCBI_spiked/GRCh38/Annotation/genes_clean.gtf
fastq_folder = /home/borzik01/liver_csc/data/fastq/GSE130473
output_directory = /home/borzik01/liver_csc/data/topout

## Assign the variable sample the list of samples to be run with tophat
samples="SRR8990162
SRR8990164
SRR8990166
SRR8990170
SRR8990172
SRR8990195
SRR8990197
SRR8990199
SRR8990201
SRR8990203
SRR8990207
SRR8990209
SRR8990210
SRR8990212"

for sample in $samples; do
echo "Running tophat on $sample..."

$TOPHAT -p 8 -o $output_directory/$sample.topout -G $GENOME $BOWTIE_INDEX $sample.fastq.gz

done

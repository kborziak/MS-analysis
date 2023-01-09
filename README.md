# NGS analysis
These scripts represent a sample of single cell RNA-seq data anaysis pipeline scripts.

The tophat_alignment.sh bash script performs a batch tophat allignment of fastq.gz files against the GRCh38 genome.

The DE analysis.R script performs differential expression analysis using EdgeR out gene count data extracted from the tophat alignments using HTSeq.

The tsne_analysis.py script performs t-SNE analysis on normalized read counts from the EdgeR analysis.

The trajectory inference.R script performs a trajectory inference analysis on normalized read counts from the EdgeR analysis.

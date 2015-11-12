# HapIso
HapIso reconstructs the haploid transcriptome of a diploid organism from long single molecule reads

Preliminary version :  the current version works per gene level

To split gene annotations by gene please run the following command :

run_splitBam.sh

To extract the list of the genes :

awk -F "gene_id" '{print $2}' genes.gtf | awk -F ";" '{print $1}' | sed 's/\"//g'  | sort | uniq >genes.LIST


To run the HapIso algorith please use hapiso_v5.py
It requires the bam file corresponding to a gene, gene coordinates and chromosome number





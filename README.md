#HapIso reconstructs the haploid transcriptome of a diploid organism from long single molecule reads

HapIso v7

This program is distributed under GPL v2. If you have questions about the license of this program, please see LICENSE in the repository.

To split bam by gene please use the following script:
run_splitBam.sh

Coordinates of the genes are stores in gene_coordinates.txt

To run the HapIso algorithm please use hapiso.py It requires the bam file corresponding to a gene, gene coordinates and chromosome number.
```
    [1] - The bam file to run HapIso
    [2] - Left boundary of the Gene
    [3] - Right boundary of the Gene
    [4] - Output file
    [5] - Chromosome
    [6] - Reference Genome in Fasta format 
```
OUTPUTS:
```
For the {input.bam} file,
the ASE will be recorded in the {input}_ASE.txt file 
the general output file will be recorded in the {input}_result.txt file 
the error corrected reads will be recorded in {input}_corrected_read.txt file. 
```








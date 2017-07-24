# HapIso reconstructs the haploid transcriptome of a diploid organism from long single molecule reads

## HapIso v.7

This program is distributed under GPL v2. If you have questions about the license of this program, please see LICENSE in the repository.

To split bam by gene please use the following script:
run_splitBam.sh

Coordinates of the genes are stores in gene_coordinates.txt

To run the HapIso algorithm please use hapiso.py It requires the bam file corresponding to a gene, gene coordinates and chromosome number
HapIso requires 6 positional arguments.

```
    [1] - The bam file to run HapIso - currently only supports one gene
    [2] - Left boundary of the Gene - the start coordniate of the gene used to extract the bam file 
    [3] - Right boundary of the Gene - the end coordniate of the gene used to extract the bam file
    [4] - Output file
    [5] - Chromosome - the chromosome of the gene used to extract the bam file
    [6] - Reference Genome in Fasta format - hg19_genome.fasta files can be found on 
```

## OUTPUTS:

```
For the {input}.bam file,
the ASE will be recorded in the {input}_ASE.txt file 
the general output file will be recorded in the {input}_result.txt file 
the error corrected reads will be recorded in {input}_corrected_read.txt file. 
```

# Contact

This software was developed by Serghei Mangul and Harry Yang. If you have any questions, please send an email to Serghei Mangul (smangul@ucla.edu) and Harry Yang (harry2416@gmail.com) . To report a bug please use github.






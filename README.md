# AlignTools

This repository contains software to generate statistics from multiple sequence alignments.

Compilation: The code can be compiled with the command: make align
Note: The Makefile may need to be edited for your specific system

To run the code, use the executable

./run_align --ali_file <alignment file name>

Outputs:

The code generates the following files:

Alignment_consensus.fa: A file giving the consensus nucleotide at each position of the alignment in fasta format.

Variant_positions.out: A list of sites in the genome at which there are genetic variants.  The format is, for example:

534 G T

here indicates that at position 534 in the genome the consensus allele is G and the second most common allele is T.

Variant_frequencies.out: The frequency of the second most common allele at each of the variant sites, expressed as a value between 0 and 1.

Variant_correlations.out: A matrix of correlations.  The position i,j in the matrix has the correlation coefficient between the occurrence of variant alleles at sites i and j.  The value 1 here would indicate that there are variants either at both sites, or at neither site.

Options:

--get_correlations 0 : Don't produce the Variant_correlations.out file

--get_frequencies 0 : Don't produce either the Variant_frequencies.out or the Variant_correlations.out files.

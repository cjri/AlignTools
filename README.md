# AlignTools

This repository contains software to generate statistics from multiple sequence alignments.

Compilation: The code can be compiled with the command: make align
Note: The Makefile may need to be edited for your specific system

To run the code, use the executable

./run_align <Method> --ali_file <alignment file name> <Options>

The method can be one of the following:

DistanceMatrix: Creates a matrix of sequence-to-sequence Hamming distances.  The matrix is output to the file Pairwise_distances.out

FilterPdiff: Filters the alignment based upon a % identity metric.

A sequence is selected from the alignment at random.  Following this, every sequence that differs at fewer than a fraction q of sites is removed, before a new sequence is selected.  The process repeats until no two remaing sequences are within the identity cutoff.

Options:

--q_cut q [Default 0.01]: Sequences with identity more than 1-q are removed [Default 99%].

Output is printed to the file Output_alignment.fa.

FilterSiteQ: Filters the alignment based upon occupancy of sites.  Used, for example, to remove sites describing rare insertion events in large alignments.

The occupancy of each site in the alignment is measured.  Sites which report a nucleotide or amino acid at less than a fraction q of the maximum are remove.

Options:

--q_cut q [Default 0.01]: Positions in the alignment with data for fewer than a fraction q of sequences are removed.

Consensus: Prints the consensus sequence of the alignment.  Output is to the file Alignment_consensus.fa.

Diversity: Calculates the diversity statistic \pi for the alignment.  Outputs to command line.

TimeSplit: Splits the alignment according to a provided list of sequence times.

Input: This routine requires a file Times.in which contains times, specified by integers, one per line, for each sequence in the alignment.  Times could be years, months, or days, for example.

Produces a series of alignment files, each containing sequences which correspond to a particular time of sequence collection.

Outputs: Where the input file is AlignmentName.xxx, sequences from each time T are output to AlignmentName_T.xxx

TimedFreqs: Calculates variant frequencies across time within an alignment.

Inupt: Requires the file Times.in, described above.

The alignment is split into sequences representing each of the times specified by the input file.  


[Work on README in progress.]




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

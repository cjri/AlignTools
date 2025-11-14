# AlignTools

This repository contains software to generate statistics from multiple sequence alignments.

**Compilation: **

The code can be compiled with the command: make align
Note that the Makefile may need to be edited for your specific system

This code uses the GSL and Eigen libraries.  On a Mac these can easily be installed using the homebrew package.

To run the code, use the executable

./run_align <Method> --ali_file <alignment file name> <Options>

The method can be one of the following:

**DistanceMatrix:** Creates a matrix of sequence-to-sequence Hamming distances.  The matrix is output to the file Pairwise_distances.out

**FilterPdiff:** Filters the alignment based upon a % identity metric.

A sequence is selected from the alignment at random.  Following this, every sequence that differs at fewer than a fraction q of sites is removed, before a new sequence is selected.  The process repeats until no two remaing sequences are within the identity cutoff.

_Options:_

--q_cut q [Default 0.01]: Sequences with identity more than 1-q are removed [Default 99%].

Output is printed to the file Output_alignment.fa.

**FilterSiteQ:** Filters the alignment based upon occupancy of sites.  Used, for example, to remove sites describing rare insertion events in large alignments.

The occupancy of each site in the alignment is measured.  Sites which report a nucleotide or amino acid at less than a fraction q of the maximum are remove.

_Options:_

--q_cut q [Default 0.01]: Positions in the alignment with data for fewer than a fraction q of sequences are removed.

**Consensus:** Prints the consensus sequence of the alignment.  Output is to the file Alignment_consensus.fa.

**Diversity:** Calculates the diversity statistic \pi for the alignment.  Outputs to command line.

**TimeSplit:** Splits the alignment according to a provided list of sequence times.

_Input:_ This routine requires a file Times.in which contains times, specified by integers, one per line, for each sequence in the alignment.  Times could be years, months, or days, for example.

Produces a series of alignment files, each containing sequences which correspond to a particular time of sequence collection.

_Outputs:_ Where the input file is AlignmentName.xxx, sequences from each time T are output to AlignmentName_T.xxx

**StringSplit:** Splits the alignment according to a provided list of sequence characteristics

_Input:_ This routine requires a file Strings.in, which contains descriptions, such as host country names, one per line, for each sequence in the alignment.

Produces a series of alignment files, each containing sequences which correspond to a particular descriptor.

_Outputs:_ Where the input file is AlignmentName.xxx, sequences corresponding to the descriptor S are output to AlignmentName_S.xxx


**TimedFreqs:** Calculates variant frequencies across time within an alignment.

_Input:_ Requires the file Times.in, described above.

The alignment is split into sequences representing each of the times specified by the input file Times.in.  Variant positions are identified in the alignment as a whole.  Next, single-nucleotide variant frequencies are calculated for each of the time-specific alignments.  Variants are derived relative to the consensus of the aligned sequences.

Variant frequency information is then output.  We identify variants that appear at a frequency of at least q_cut in a time-specific alignment containing at least n_cut sequences.  

_Output:_ Variant counts are then output to the file Variant_trajectories.out, which has the format:

Position in genome	Consensus	Alternative	#Times	(Time	#A	#C	#G	#T	Total)

for nucleotide sequences; for protein sequences the number of each amino acid is shown.

_Options:_

--q_cut q [Default 0.01]: Minimum frequency at which a variant can be identified and reported.

--n_cut n [Default 10]: Minimum number of sequences required to identify a variant in a specific alignment.

--n_reps r [Default 1]: Minumum number of occasions on which a variant needs to be observed before it is reported.

**Random:** Generates random sequences matching the properties of the alignment

We generate variant frequencies, including single nucleotide polymorphisms, and deletions (insertions are presupposed to have been incorporated into the alignment), across the sequences in the alignment.  We then generate a correlation matrix, describing pairwise correlations between these variants.  Random sequences are then generated in a manner that reflects the variant frequencies and the correlations between them.

_Output:_ Generated sequences are reported either as a series of single-nucleotide variant positions or as complete sequences.  

_Options:_

--generate <number> : Number of random sequences to generate.

--output <type> [Default: Sparse] : Output format for random sequences:

Sparse: Outputs comprise nucleotide positions at which single nucleotide variants arise.  Deletions are not explicitly reported, though the locations of variants that fall in deletions are not reported.  This format was designed for sparse storage of potentially large numbers of sequences, generated from alignments of viral genomes, where indels are often rare.

FASTA: Outputs comprise sequences in FASTA format.  Deletions are explicitly included in sequences.  The locations of variants is not reported, but could be derived from a comparison of sequences if needed.
	



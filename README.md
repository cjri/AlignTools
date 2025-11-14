# AlignTools

This repository contains software to generate statistics from multiple sequence alignments.

Compilation: The code can be compiled with the command: make align
Note: The Makefile may need to be edited for your specific system

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


**TimedFreqs:** Calculates variant frequencies across time within an alignment.

_Input:_ Requires the file Times.in, described above.

The alignment is split into sequences representing each of the times specified by the input file Times.in.  Variant positions are identified in the alignment as a whole.  Next, variant frequencies are calculated for each of the time-specific alignments.

Variant frequency information is then output.  We identify variants that appear at a frequency of at least q_cut in a time-specific alignment containing at least n_cut sequences.  

_Output:_ Variant counts are then output to the file Variant_trajectories.out, which has the format:

Position in genome	Consensus	Alternative	#Times	(Time	#A	#C	#G	#T	Total)

for nucleotide sequences; for protein sequences the number of each amino acid is shown.

_Options:_

--q_cut q [Default 0.01]: Minimum frequency at which a variant can be identified and reported.

--n_cut n [Default 10]: Minimum number of sequences required to identify a variant in a specific alignment.


#include "aligntools.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void MakeDistanceMatrix (run_params& p, int output, vector<string>& seqs, vector<string>& names, vector< vector<int> >& seqdists) {
    string all_consensus;
    FindConsensus(all_consensus,seqs);
    vector<sparseseq> variants;
    FindSVariants (variants,all_consensus,seqs);
    //Thought: Do we want to check for ambiguous nucleotides.  Currently no.
    
    FindPairwiseDistances (variants,seqs,seqdists);
    
    if (output==1) {
        ofstream pair_file;
        pair_file.open("Pairwise_distances.out");
        for (int i=0;i<seqdists.size();i++) {
            for (int j=0;j<seqdists[i].size();j++) {
                pair_file << seqdists[i][j] << " ";
            }
            pair_file << "\n";
            if (p.dist_cut>0) {
                vector< vector<int> > subsets;
                GetSubsetsIJ (p,names,seqdists,subsets);
            }
        }
    }
}

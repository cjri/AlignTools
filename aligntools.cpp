//Program to calculate predicted sigma values for loci in a model system, where data is generated by a sampling process
//This code uses a deterministic method to run the underlying population distribution, and samples it at low frequency

#include <iostream>
#include <sstream>
#include <vector>
#include <list>
#include <deque>
using namespace std;

#include "aligntools.h"
#include "utilities.h"

int main(int argc, const char **argv){

    //Code to read in an alignment.
    //Identifies variants and will make covariance matrix between these sites.
    
	run_params p;
	GetParameters(p,argc,argv);

    if (p.dismat==1) {
        //Create distance matrix between sequences
        vector<string> seqs;
        vector<string> names;
        ReadFastaAli(p,seqs,names);
        CheckBaseCase(seqs);
        string all_consensus;
        FindConsensus(all_consensus,seqs);
        vector<sparseseq> variants;
        FindSVariants (variants,all_consensus,seqs);
        //Thought: Do we want to check for ambiguous nucleotides.  Currently no.
        
        vector< vector<int> > seqdists;
        FindPairwiseDistances (variants,seqs,seqdists);
        for (int i=0;i<seqdists.size();i++) {
            
            for (int j=0;j<seqdists[i].size();j++) {
                cout << seqdists[i][j] << " ";
            }
            cout << "\n";
        }
        
        if (p.dist_cut>0) {
            cout << "Here\n";
            vector< vector<int> > subsets;
            GetSubsetsIJ (p,names,seqdists,subsets);
        }
    }
    
	//Open alignment
    vector<site> ali_stats;
    ReadVariants (p,ali_stats);
            
    //Next step - find variants
    vector<int> var_positions;
    FindVariants (ali_stats,var_positions);

    //Consensus sequence
    vector<string> consensus;
    GetConsensus(ali_stats,consensus);
    
    vector<string> second=consensus;
    CalculateFrequencies (ali_stats,second);
    
    ofstream vars_file;
    vars_file.open("Variant_positions.out");
    for (int i=0;i<var_positions.size();i++) {
        vars_file << var_positions[i]+1 << " " << consensus[var_positions[i]] << " " << second[var_positions[i]] << "\n";
    }
    
    //Get frequencies - binary
    if (p.get_frequencies==1) {
        ofstream freqs_file;
        freqs_file.open("Variant_frequencies.out");
        for (int i=0;i<var_positions.size();i++) {
            freqs_file << ali_stats[var_positions[i]].freq << "\n";
        }
    }
    
    //Make matrices - might need to read again?  Or do from the whole alignment.
    
    if (p.get_correlations==1) {
        vector<pr> pairs;
        MakeInitialPairs (var_positions,pairs);
        ConstructPairs (p,second,pairs);

        //Find correlations
        FindCorrelations (ali_stats,pairs);
        ofstream correl_file;
        correl_file.open("Variant_correlations.out");
        int index=0;
        for (int i=0;i<var_positions.size();i++) {
            for (int j=0;j<var_positions.size();j++) {
                correl_file << pairs[index].correl << " ";
                index++;
            }
            correl_file << "\n";
        }
    }
    
    
    
    return 0;
}
	
	

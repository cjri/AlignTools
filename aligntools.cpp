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
	
	gsl_rng_env_setup();
	gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
	gsl_rng_set (rgen, p.seed);

	//Open alignment
    vector<site> ali_stats;
    ReadVariants (p,ali_stats);
            
    //Next step - find variants
    vector<int> var_positions;
    FindVariants (ali_stats,var_positions);

    //Consensus sequence
    vector<string> consensus;
    GetConsensus(ali_stats,consensus);
    
    //Get frequencies - binary
    vector<string> second=consensus;
    CalculateFrequencies (ali_stats,second);
    
    //Make matrices - might need to read again?  Or do from the whole alignment.
    
    vector<pr> pairs;
    MakeInitialPairs (var_positions,pairs);
    ConstructPairs (p,second,pairs);

    //Find correlations
    FindCorrelations (ali_stats,pairs);
    
    ofstream vars_file;
    vars_file.open("Variant_positions.out");
    for (int i=0;i<var_positions.size();i++) {
        vars_file << var_positions[i]+1 << " " << consensus[var_positions[i]] << " " << second[var_positions[i]] << "\n";
    }
    
    ofstream freqs_file;
    freqs_file.open("Variant_frequencies.out");
    for (int i=0;i<var_positions.size();i++) {
        freqs_file << ali_stats[var_positions[i]].freq << "\n";
    }
    
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
    
    return 0;
}
	
	

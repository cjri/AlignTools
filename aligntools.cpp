#include <iostream>
#include <sstream>
#include <vector>
#include <list>
using namespace std;

#include "aligntools.h"
#include "distance_matrix.h"
#include "diversity.h"
#include "editing.h"
#include "io.h"
#include "rgen.h"
#include "stringsplit.h"
#include "timesplit.h"
#include "utilities.h"
#include <Eigen/Dense>

int main(int argc, const char **argv){

    //Code to read in an alignment.
    //Identifies variants and will make covariance matrix between these sites.
    
    //Initialise random number generator
    
	run_params p;
	GetParameters(p,argc,argv);
    p.seed=(int) time(NULL);
    gsl_rng_env_setup();
    gsl_rng *rgen = gsl_rng_alloc (gsl_rng_taus);
    gsl_rng_set (rgen, p.seed);
    
    vector<string> seqs;
    vector<string> names;
    ReadFastaAli(p,seqs,names);
    CheckBaseCase(seqs);
    if (p.error==1) {
        cout << "Terminate with error\n";
        return 0;
    }
    vector<char> alphabet;
    GetAlignmentType(p,alphabet,seqs);
    
    //To add: Replace N in a nucleotide alignment with - if this is being used for insertions?
    
    
    
    if (p.method.compare("DistanceMatrix")==0) {
        //Updated
        vector< vector<int> > seqdists;
        MakeDistanceMatrix (p,1,alphabet,seqs,names,seqdists);
        return 0;
    }
    
    if (p.method.compare("FilterPdiff")==0) {
        //Updated
        FilterPDiff (p,seqs,names,alphabet,rgen);
        return 0;
    }

    //Get alignment statistics
    //vector<site> ali_stats;
    vector<site2> ali_stats2;
    GetAliStats2 (p,seqs,alphabet,ali_stats2);
    
    if (p.verb==1) {
        cout << "Alignment statistics\n";
        for (int i=0;i<ali_stats2.size();i++) {
            for (int j=0;j<alphabet.size();j++) {
                cout << alphabet[j] << " " << ali_stats2[i].counts[j] << " ";
            }
            cout << ali_stats2[i].N << "\n";
        }
    }
    
    if (p.method.compare("FilterSiteQ")==0) {
        //Updated
        FilterAlignmentQ (p,ali_stats2,seqs,names);
        return 0;
    }

    
    //Consensus sequence
    vector<string> consensus;
    GetConsensus2(ali_stats2,alphabet,consensus);
    int seq_length=consensus.size();
    
    if (p.method.compare("Consensus")==0) {
        //Updated
        return 0;
    }

    //Next step - variant positions
    vector<int> var_positions;
    FindVariants (ali_stats2,var_positions);
    
    //cout << "Positions " << var_positions.size() << " " << consensus.size() << "\n";
    
    if (p.method.compare("Diversity")==0) {
        //Updated
        GetPiDiversity (p,seq_length,alphabet,seqs,names);
        return 0;
    }

    if (p.method.compare("TimeSplit")==0) {
        TDSplit (p,consensus,var_positions,ali_stats2,names,seqs);
        return 0;
    }
    
    if (p.method.compare("StringSplit")==0) {
        SSplit (p,consensus,var_positions,ali_stats2,names,seqs);
        return 0;
    }

    
    if (p.method.compare("Random")==0) {
        //To update
        GenerateRandomSequences (p,seq_length,consensus,alphabet,var_positions,ali_stats2,seqs,rgen);
        return 0;
    }
        
    if (p.method.compare("TimedFreqs")==0) {
        //Updated
        if (p.type==1) {
            GetTDNucleotideCounts (p,consensus,alphabet,var_positions,ali_stats2,seqs);
        }
        return 0;
        
    } else {
        cout << "Instructions:\n";
        cout << "./run_align DistanceMatrix <flags> : Calculates matrix of distances between sequences.\n";
        cout << "\n";
        cout << "./run_align Diversity <flags> : Calculates pi diversity for sequences in the alignment.\n";
        cout << "\n";
        cout << "./run_align Random <flags> : Generates random sequences.\n";
        cout << "Flags are:\n";
        cout << " --generate <number> : Number of random sequences to generate\n";
        cout << " --verb <number> : Write out variant details to file\n";
        cout << " --output <type> : Output format for random files:\n";
        cout << "       Sparse: [Default] Output positions of variants w.r.t. consensus\n";
        cout << "       FASTA: Output as .fasta file format\n";
        cout << "       Binary: Output as binary string at variant positions\n";
        cout << "\n";
        cout << "./run_align FilterSiteQ <flags> : Removes positions in the genome which don't contain {A, C, G, T} nucleotides for many of the sequences.\n";
        cout << "Flags are:\n";
        cout << " --q_cut <frequency> : [Default 0.1] Require this fraction of sites to have an {A,C,G,T} nucleotide\n";
        cout << "Output to Output_alignment.fa\n";
        cout << "\n";
        cout << "./run_align TimedFreqs <flags> : Separates sequences in an alignment by time and produces a record of variant frequency over time.\n";
        cout << "Flags are:\n";
        cout << " --q_cut <frequency> : [Default 0.1] Minimum minor allele frequency to report\n";
        cout << " --n_cut <value> : [Default 10] Minimum read depth when calling variant allele\n";
        cout << " --n_reps <type> : [Default 1] Minimum number of times variant observed at given frequency and read depth:\n";
        cout << "\n";
        cout << "./run_align TimeSplit <flags> : Splits alignment into time-specific sub-alignments and printes these alignments.\n";
        cout << "Flags are:\n";
        cout << " --ali_file <file> : Multiple sequence alignment file in .fasta format\n";
        cout << "\n";

    }
    
    return 0;
}
	
	

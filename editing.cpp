#include "aligntools.h"
#include "editing.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void FilterAlignmentQ (run_params& p, vector<site>& ali_stats, vector<string>& seqs, vector<string>& names) {
    //Filter the genome sequence alignment keeping only positions at which at least a fraction p.qq_cut of the sequences have an {A,C,G,T} nucleotide.
    int max=0;
    for (int i=0;i<ali_stats.size();i++) {
        if (ali_stats[i].N>max) {
            max=ali_stats[i].N;
        }
    }
    //cout << "Max depth " << max << "\n";
    int threshold=max*p.qq_cut;
    //cout << "Threshold " << threshold << "\n";
    vector<int> keep;
    for (int i=0;i<ali_stats.size();i++) {
        if (ali_stats[i].A+ali_stats[i].C+ali_stats[i].G+ali_stats[i].T>=threshold) {
            keep.push_back(i);
            //cout << i << " " << ali_stats[i].A+ali_stats[i].C+ali_stats[i].G+ali_stats[i].T << "\n";
        }
    }
    //cout << "Sites to keep " << keep.size() << " of " << ali_stats.size() << "\n";
    OutputAlignmentFiltered (names,seqs,keep);
}

void FilterPDiff (run_params& p, vector<string>& seqs, vector<string>& names, gsl_rng *rgen) {
    //Filter the genome sequence alignment, repeatedly choosing a random sequence, then removing anything with more than p.qq_cut fraction difference
    //Makes a smaller alignment covering the overall diversity of the sequences
    string all_consensus;
    FindConsensus(all_consensus,seqs);
    vector<sparseseq> variants;
    FindSVariants (variants,all_consensus,seqs);
    //cout << variants.size() << "\n";
    //cout << "Length " << seqs[0].length() << "\n";
    
    p.qq_cut=0.01;
    int threshold=seqs[0].length()*p.qq_cut;
    //Keep track of sequences
    vector<int> done;
    vector<int> index;
    for (int i=0;i<seqs.size();i++) {
        done.push_back(0);
    }
    
    vector<int> sample;
    int donezero=1;
    while (donezero==1) {
        index.clear();
        for (int i=0;i<done.size();i++) {
            if (done[i]==0) {
                index.push_back(i);
            }
        }
        //cout << "Index size " << index.size() << "\n";
        int selected=floor(gsl_rng_uniform(rgen)*index.size()+0.5)-1;
        selected=index[selected];
        done[selected]=1;
        sample.push_back(selected);
        
        //Find everything within threshold and remove
        for (int i=0;i<variants.size();i++) {
            if (done[i]==0) {
                //Find unique difference positions between selected and i
                vector<int> uniq;
                for (int k=0;k<variants[i].locus.size();k++) {
                    uniq.push_back(variants[i].locus[k]);
                }
                for (int k=0;k<variants[selected].locus.size();k++) {
                    uniq.push_back(variants[selected].locus[k]);
                }
                sort(uniq.begin(),uniq.end());
                uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
                
                //Run comparison over uniq positions
                int dist=0;
                int kk=0;
                while (kk<uniq.size()&&dist<threshold) {
                    if (seqs[i][uniq[kk]]=='A'||seqs[i][uniq[kk]]=='C'||seqs[i][uniq[kk]]=='G'||seqs[i][uniq[kk]]=='T') {
                        if (seqs[selected][uniq[kk]]=='A'||seqs[selected][uniq[kk]]=='C'||seqs[selected][uniq[kk]]=='G'||seqs[selected][uniq[kk]]=='T') {
                            if (seqs[i][uniq[kk]]!=seqs[selected][uniq[kk]]) {
                                dist++;
                            }
                        }
                    }
                    kk++;
                }
                if (dist<threshold) {
                    done[i]=1;
                }
            }
        }
        donezero=0;
        for (int i=0;i<done.size();i++) {
            if (done[i]==0) {
                donezero=1;
                break;
            }
        }
        
        
    }
    sort(sample.begin(),sample.end());
    sample.erase(unique(sample.begin(),sample.end()),sample.end());

    /*cout << "Sample\n";
    for (int i=0;i<sample.size();i++) {
        cout << sample[i] << " ";
    }
    cout << "\n";*/
    OutputAlignmentSFiltered (names,seqs,sample);
    
    

}

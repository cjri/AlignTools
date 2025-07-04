#include "aligntools.h"
#include "io.h"
#include "rgen.h"
#include <iostream>
#include <string>
#include <cstring>
#include <random>

void ReadFastaAli (run_params p, vector<string>& seqs, vector<string>& names) {
    ifstream ali_file;
    ali_file.open(p.ali_file.c_str());
    string seq;
    string str;
    for (int i=0;i<1000000;i++) {
        if (!(ali_file >> str)) break;
        if (str.at(0)=='>') {
            str.erase(0,1);
            names.push_back(str);
            if (seq.size()>0) {
                seqs.push_back(seq);
                seq.clear();
            }
        } else {
            seq=seq+str;
        }
    }
    if (seq.size()>0) {
        seqs.push_back(seq);
    }
}

void ReadTimes (vector<int>& times, vector<int>& times_uniq) {
    ifstream t_file;
    t_file.open("Times.in");
    int t;
    for (int i=0;i<1000000;i++) {
        if (!(t_file >> t)) break;
        times.push_back(t);
    }
    times_uniq=times;
    sort(times_uniq.begin(),times_uniq.end());
    times_uniq.erase(unique(times_uniq.begin(),times_uniq.end()),times_uniq.end());
}

void PrintVariantPositions (vector<int>& var_positions, vector<string>& consensus, vector<string>& second) {
    ofstream vars_file;
    vars_file.open("Variant_positions.out");
    for (int i=0;i<var_positions.size();i++) {
        vars_file << var_positions[i]+1 << " " << consensus[var_positions[i]] << " " << second[var_positions[i]] << "\n";
    }
}

void PrintFrequencies (vector<double>& frequencies) {
    ofstream freqs_file;
    freqs_file.open("Variant_frequencies.out");
    for (int i=0;i<frequencies.size();i++) {
        freqs_file << frequencies[i] << "\n";
    }
}

void PrintCorrelations (vector< vector<double> >& correls) {
    ofstream correl_file;
    correl_file.open("Variant_correlations.out");
    for (int i=0;i<correls.size();i++) {
        vector<double> c;
        for (int j=0;j<correls[i].size();j++) {
            correl_file << correls[i][j] << " ";
        }
        correl_file << "\n";
    }
}

void OutputBitstrings (run_params& p, vector< vector<int> >& bitstrings, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<delet> deletions) {
    if (p.output.compare("Sparse")==0) {
        OutputBitstringsSparse(var_positions,deletions,bitstrings);
    }
    if (p.output.compare("FASTA")==0) {
        OutputBitstringsFasta(var_positions,consensus,second,deletions,bitstrings);
    }
    if (p.output.compare("Binary")==0) {
        OutputBitstringsBinary(var_positions,deletions,bitstrings);
    }
    
}

void OutputBitstringsBinary (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings) {
    //Replace deletions by zeros
    for (int d=0;d<deletions.size();d++) {
        for (int i=0;i<bitstrings.size();i++) {
            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(0.0, 1.0);
            double r = dis(gen);
            if (r<deletions[i].freq) {
                for (int j=0;j<var_positions.size();j++) {
                    for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                        if (var_positions[j]==k) {
                            bitstrings[i][j]=0;
                        }
                    }
                }
            }
        }
    }
    for (int i=0;i<bitstrings.size();i++) {
        cout << "Sample " << i+1 << " ";
        for (int j=0;j<bitstrings[i].size();j++) {
            cout << bitstrings[i][j] << " ";
        }
        cout << "\n";
    }
}

void OutputBitstringsSparse (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings) {
    //Sparse variant format
    //Replace deletions by zeros
    for (int d=0;d<deletions.size();d++) {
        for (int i=0;i<bitstrings.size();i++) {
            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(0.0, 1.0);
            double r = dis(gen);
            if (r<deletions[i].freq) {
                for (int j=0;j<var_positions.size();j++) {
                    for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                        if (var_positions[j]==k) {
                            bitstrings[i][j]=0;
                        }
                    }
                }
            }
        }
    }
    for (int i=0;i<bitstrings.size();i++) {
        cout << "Sample " << i+1 << " ";
        for (int j=0;j<bitstrings[i].size();j++) {
            if (bitstrings[i][j]==1) {
                cout << var_positions[j] << " ";
            }
        }
        cout << "\n";
    }
}

void OutputBitstringsFasta (const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, const vector<delet>& deletions, vector< vector<int> >& bitstrings) {
    //Replace deletions by negative 1: Flag for a '-'
    for (int d=0;d<deletions.size();d++) {
        for (int i=0;i<bitstrings.size();i++) {
            random_device rd;
            mt19937 gen(rd());
            uniform_real_distribution<> dis(0.0, 1.0);
            double r = dis(gen);
            if (r<deletions[i].freq) {
                for (int j=0;j<var_positions.size();j++) {
                    for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                        if (var_positions[j]==k) {
                            bitstrings[i][j]=-1;
                        }
                    }
                }
            }
        }
    }

    for (int i=0;i<bitstrings.size();i++) {
        int index=0;
        cout << ">Sample" << i+1 << "\n";
        for (int j=0;j<bitstrings[i].size();j++) {
            while (index<var_positions[j]) {
                cout << consensus[index];
                index++;
            }
            if (bitstrings[i][j]==0) {
                cout << consensus[index];
            } else if (bitstrings[i][j]==-1){
                cout << "-";
            } else {
                cout << second[index];
            }
            index++;
        }
        while (index<consensus.size()) {
            cout << consensus[index];
            index++;
        }
        cout<< "\n";
    }
}

void OutputNucleotideCountsTime (run_params& p, vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site> >& ali_stats_t)  {
    ofstream nct_file;
    nct_file.open("Variant_trajectories.out");
    for (int i=0;i<var_positions.size();i++) { //All variants from all times
        //Check for maximum variant frequency
        int include=0;
        for (int t=0;t<times_uniq.size();t++) {
            if (ali_stats_t[t][var_positions[i]].freq>=p.qq_cut&&ali_stats_t[t][var_positions[i]].N>=p.n_cut) {
                include++;
            }
        }
        if (include>=p.n_reps) {
            nct_file << var_positions[i] << " " << consensus[var_positions[i]] << " " << second[var_positions[i]] << " ";
            for (int t=0;t<times_uniq.size();t++) {
                nct_file << times_uniq[t] << " " << ali_stats_t[t][var_positions[i]].A << " " << ali_stats_t[t][var_positions[i]].C << " " << ali_stats_t[t][var_positions[i]].G << " " << ali_stats_t[t][var_positions[i]].T << " " << ali_stats_t[t][var_positions[i]].N << " ";
            }
            nct_file << "\n";
        }
    }
}

void OutputVariantFrequencies (vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site> >& ali_stats_t)  {
    ofstream ncf_file;
    ncf_file.open("Variant_frequencies.out");
    for (int i=0;i<var_positions.size();i++) { //All variants from all times
        ncf_file << var_positions[i] << " " << consensus[var_positions[i]] << " " << second[var_positions[i]] << " ";
        for (int t=0;t<times_uniq.size();t++) {
            ncf_file << times_uniq[t] << " " << ali_stats_t[t][var_positions[i]].freq << " ";
        }
        ncf_file << "\n";
    }
}

#include "aligntools.h"
#include "io.h"
#include "rgen.h"
#include <iostream>
#include <string>
#include <cstring>
#include <random>

void ReadFastaAli (run_params& p, vector<string>& seqs, vector<string>& names) {
    ifstream ali_file;
    ali_file.open(p.ali_file.c_str());
    string seq;
    string str;
    while (getline(ali_file, str)) {
        if (str.empty()) continue;
        if (str.at(0)=='>') {
            str.erase(0,1);
            names.push_back(str);
            if (!seq.empty()) {
                seqs.push_back(seq);
                seq.clear();
            }
        } else {
            seq=seq+str;

        }
    }
    if (seqs.size()==0) {
        p.error=1;
        cout << "Error: Alignment file not found.  Specify using --ali_file <filename>\n";
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

void OutputBitstrings (run_params& p, vector< vector<int> >& bitstrings, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<delet> deletions, vector< vector<int> >& denovo) {
    if (p.output.compare("Sparse")==0) {
        OutputBitstringsSparse(var_positions,deletions,bitstrings,denovo);
    }
    if (p.output.compare("FASTA")==0) {
        OutputBitstringsFasta(var_positions,consensus,second,deletions,bitstrings,denovo);
    }
    if (p.output.compare("Binary")==0) {
        OutputBitstringsBinary(var_positions,deletions,bitstrings,denovo);
    }
    
}

void OutputBitstringsBinary (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo) {
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
                //Sort out de novo mutations
                for (int j=0;j<denovo[i].size();j++) {
                    for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                        if (denovo[i][j]==k) {
                            denovo[i][j]=-1;
                        }
                    }
                }
            }
        }
    }
    EditDenovo(denovo);
    //Insert de novo mutations
    vector<int> all_denovo;
    GetDenovoSites (denovo,all_denovo);
    
    OutputSites (var_positions,all_denovo,bitstrings);
    
    for (int i=0;i<bitstrings.size();i++) {
        int indexb=0;
        int indexd=0;
        while (indexb<bitstrings.size()||indexd<all_denovo.size()) {
            if (var_positions[indexb]<all_denovo[indexd]) {
                cout << bitstrings[i][indexb] << " ";
                indexb++;
            } else {
                int found=0;
                for (int j=0;j<denovo[i].size();j++) {
                    if (denovo[i][j]==all_denovo[indexd]) {
                        found=1;
                        break;
                    }
                }
                if (found==1) {
                    cout << "1 ";
                } else {
                    cout << "0 ";
                }
                indexd++;
            }
        }
        cout << "\n";
    }
}

void GetDenovoSites (vector< vector<int> >& denovo, vector<int>& all_denovo) {
    for (int i=0;i<denovo.size();i++) {
        for(int j=0;j<denovo[i].size();j++) {
            all_denovo.push_back(denovo[i][j]);
        }
    }
    if (all_denovo.size()>0) {
        sort(all_denovo.begin(),all_denovo.end());
        all_denovo.erase(unique(all_denovo.begin(),all_denovo.end()),all_denovo.end());
    }
}

void OutputSites (const vector<int>& var_positions, vector<int>& all_denovo, vector< vector<int> >& bitstrings) {
    cout << "Sites\n";
    int indexb=0;
    int indexd=0;
    while (indexb<bitstrings.size()||indexd<all_denovo.size()) {
        if (var_positions[indexb]<all_denovo[indexd]) {
            cout << var_positions[indexb] << " ";
            indexb++;
        } else {
            cout << all_denovo[indexd] << " ";
            indexd++;
        }
    }
    cout << "\n";
}

void OutputBitstringsSparse (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo) {
    cout << "Sparse output\n";
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
                //Sort out de novo mutations
                for (int j=0;j<denovo[i].size();j++) {
                    for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                        if (denovo[i][j]==k) {
                            cout << "Delete " << denovo[i][j] << "\n";
                            denovo[i][j]=-1;
                        }
                    }
                }
            }
        }
    }
    EditDenovo(denovo);

    //Insert de novo mutations
    vector<int> all_denovo;
    GetDenovoSites (denovo,all_denovo);

    cout << "Denovo sites " << all_denovo.size() << "\n";

    
    for (int i=0;i<bitstrings.size();i++) {
        cout << "Sample " << i+1 << " ";
        int indexb=0;
        int indexd=0;
        while (indexb<var_positions.size()||indexd<all_denovo.size()) {
            //Find the next value
            int min=100000;
            int fromb=0;
            if (indexb<var_positions.size()&&var_positions[indexb]<min) {
                min=var_positions[indexb];
                fromb=1;
            }
            if (indexd<all_denovo.size()&&all_denovo[indexd]<min) {
                min=all_denovo[indexd];
                fromb=0;
            }

            if (fromb==1) {
                if (bitstrings[i][indexb]==1) {
                    cout << var_positions[indexb] << " ";
                }
                indexb++;
            } else {
                for (int j=0;j<denovo[i].size();j++) {
                    if (denovo[i][j]==all_denovo[indexd]) {
                        cout << all_denovo[indexd] << " ";
                        break;
                    }
                }
                indexd++;
            }
        }
        cout << "\n";
    }
}

void OutputBitstringsFasta (const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo) {

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
            //Do the same for de novo mutations in this sequence
            for (int j=0;j<denovo[i].size();j++) {
                for (int k=deletions[d].start;k<deletions[d].start+deletions[d].length;k++) {
                    if (denovo[i][j]==k) {
                        denovo[i][j]=-1;
                    }
                }
            }
        }
    }

    //Remove negative denovo mutations
    EditDenovo(denovo);
    
    //Insert de novo mutations
    vector<int> all_denovo;
    GetDenovoSites (denovo,all_denovo);

    cout << "All denovo size " << all_denovo.size() << "\n";
    
    cout << "Var_positions size " << var_positions.size() << "\n";
    
    
    
    for (int i=0;i<bitstrings.size();i++) {
        int indexb=0;
        int indexd=0;
        int indexs=0;
        cout << ">Sample" << i+1 << "\n";
        while (indexb<var_positions.size()||indexd<all_denovo.size()) {
            int min=100000;
            int fromb=0;
            if (indexb<var_positions.size()&&var_positions[indexb]<min) {
                min=var_positions[indexb];
                fromb=1;
            }
            if (indexd<all_denovo.size()&&all_denovo[indexd]<min) {
                min=all_denovo[indexd];
                fromb=0;
            }
            while (indexs<min) {
                cout << consensus[indexs];
                indexs++;
            }
            
            if (fromb==1) {
                if (bitstrings[i][indexb]==0) {
                    cout << consensus[indexs];
                } else if (bitstrings[i][indexb]==-1) {
                    cout << "-";
                } else {
                    cout << second[indexs];
                }
                indexs++;
                indexb++;
            } else {
                int found=0;
                for (int j=0;j<denovo[i].size();j++) {
                    if (denovo[i][j]==all_denovo[indexd]) {
                        found=1;
                        break;
                    }
                }
                if (found==0) {
                    cout << consensus[indexs];
                } else {
                    //New denovo mutation
                    if (consensus[indexs]=="A") {
                        cout << "C";
                    }
                    if (consensus[indexs]=="C") {
                        cout << "G";
                    }
                    if (consensus[indexs]=="G") {
                        cout << "T";
                    }
                    if (consensus[indexs]=="T") {
                        cout << "A";
                    }
                    
                }
                indexs++;
                indexd++;
            }
        }
        cout << "\n";
    }
    
}

void EditDenovo (vector< vector<int> >& denovo) {
    vector< vector<int> > dn_new;
    for (int i=0;i<denovo.size();i++) {
        vector<int> replace;
        for (int j=0;j<denovo[i].size();j++) {
            if (denovo[i][j]>=0) {
                replace.push_back(denovo[i][j]);
            }
        }
        dn_new.push_back(replace);
    }
    denovo=dn_new;
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

void OutputAlignmentFiltered (vector<string>& names, vector<string>& seqs, vector<int>& keep) {
    //Filtered by nucleotide positions
    ofstream a_file;
    a_file.open("Output_alignment.fa");
    for (int i=0;i<names.size();i++) {
        a_file << ">" << names[i] << "\n";
        for (int j=0;j<keep.size();j++) {
            if (seqs[i][keep[j]]=='N') {
                a_file << "-";
            } else {
                a_file << seqs[i][keep[j]];
            }
        }
        a_file << "\n";
    }
}

void OutputAlignmentSFiltered (vector<string>& names, vector<string>& seqs, vector<int>& sample) {
    //Filtered by sequences
    ofstream a_file;
    a_file.open("Output_alignment.fa");
    for (int i=0;i<sample.size();i++) {
        a_file << ">" << names[sample[i]] << "\n";
        for (int j=0;j<seqs[sample[i]].size();j++) {
            if (seqs[sample[i]][j]=='N') {
                a_file << "-";
            } else {
                a_file << seqs[sample[i]][j];
            }
        }
        a_file << "\n";
    }
}

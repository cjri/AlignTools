#include "aligntools.h"
#include "io.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>
#include <unordered_map>
#include <algorithm>


void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
    p.method="Instructions";
	p.ali_file="Align.fa";
    p.dismat=0;
    p.get_positions=1;
    p.get_frequencies=1;
    p.get_correlations=1;
    p.n_generations=0;
    p.cutoff=0.0;
    p.qq_cut=0.01;
    p.n_cut=10;
    p.n_reps=1;
    p.denovo=0;
    p.output="Sparse"; //Options FASTA, Binary
    p.verb=0;
    p.error=0;
    p.type=0;
    p.method=argv[x];
    cout << "Method " << p.method << "\n";
    x++;
	while (x < argc && (argv[x][0]=='-')) {
		p_switch=argv[x];
		if (p_switch.compare("--ali_file")==0) {
			x++;
			p.ali_file=argv[x];
		} else if (p_switch.compare("--get_frequencies")==0) {
            x++;
            p.get_frequencies=atoi(argv[x]);
        } else if (p_switch.compare("--get_correlations")==0) {
            x++;
            p.get_correlations=atoi(argv[x]);
        } else if (p_switch.compare("--distances")==0) {
            x++;
            p.dismat=atoi(argv[x]);
        } else if (p_switch.compare("--dist_cut")==0) {
            x++;
            p.dist_cut=atoi(argv[x]);
        } else if (p_switch.compare("--generate")==0) {
            x++;
            p.n_generations=atoi(argv[x]);
        } else if (p_switch.compare("--verb")==0) {
            x++;
            p.verb=atoi(argv[x]);
        } else if (p_switch.compare("--q_nocorrel")==0) {
            x++;
            p.cutoff=atof(argv[x]);
        } else if (p_switch.compare("--q_cut")==0) {
            x++;
            p.qq_cut=atof(argv[x]);
        } else if (p_switch.compare("--n_cut")==0) {
            x++;
            p.n_cut=atoi(argv[x]);
        } else if (p_switch.compare("--n_reps")==0) {
            x++;
            p.n_reps=atoi(argv[x]);
        } else if (p_switch.compare("--denovo")==0) {
            x++;
            p.denovo=atoi(argv[x]);

        } else if (p_switch.compare("--output")==0) {
            x++;
            p.output=argv[x];
        } else {
			cout << "Incorrect usage\n ";
            cout << p_switch << "\n";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
    if (p.get_frequencies==0) {
        p.get_correlations=0;
    }
}


void CheckBaseCase (vector<string>& seqs) {
    for (int i=0;i<seqs.size();i++) {
        transform(seqs[i].begin(), seqs[i].end(), seqs[i].begin(),
                  ::toupper);
    }
}

void GetAlignmentType (run_params& p, vector<char>& alphabet, vector<string>& seqs) {
    //Uses first sequence.  Protein or nucleotide
    vector<char> nuc_chars = {'A','C','G','T'};
    vector<char> aa_chars = {'A','C','D','E','F','G','H','I','L','K','M','N','P','Q','R','S','T','V','W','Y'};
    double nuc_count = 0;
    //double aa_count  = 0;
    double valid_count = 0;
    for (int i=0;i<seqs[0].size();i++) {
        int to_add=0;
        char charToFind = seqs[0][i];
        auto it1 = std::find(nuc_chars.begin(), nuc_chars.end(), charToFind);
        auto it2 = std::find(aa_chars.begin(), aa_chars.end(), charToFind);
        if (it1 != nuc_chars.end()) {
            nuc_count++;
            to_add=1;
        }
        if (it2 != aa_chars.end()) {
            //aa_count++;
            to_add=1;
        }
        if (to_add==1) {
            valid_count++;
        }
    }
    if (nuc_count>(0.85*valid_count)) {
        p.type=0;
        alphabet=nuc_chars;
    } else {
        p.type=1;
        alphabet=aa_chars;
    }
    if (p.type==0) {
        cout << "Detected nucleotide alignment\n";
    }
    if (p.type==1) {
        cout << "Detected amino acid alignment\n";
    }
}

/*void FindConsensus (string& consensus, vector<string>& seqs) {
    consensus=seqs[0];
    if (seqs.size()>1) {
        int nA=0;
        int nC=0;
        int nG=0;
        int nT=0;
        for (int pos=0;pos<seqs[0].size();pos++) {
            nA=0;
            nC=0;
            nG=0;
            nT=0;
            for (int seq=0;seq<seqs.size();seq++) {
                if (seqs[seq][pos]=='A') {
                    nA++;
                }
                if (seqs[seq][pos]=='C') {
                    nC++;
                }
                if (seqs[seq][pos]=='G') {
                    nG++;
                }
                if (seqs[seq][pos]=='T') {
                    nT++;
                }
            }
            int max=nA;
            consensus[pos]='A';
            if (nC>max) {
                max=nC;
                consensus[pos]='C';
            }
            if (nG>max) {
                max=nG;
                consensus[pos]='G';
            }
            if (nT>max) {
                consensus[pos]='T';
                max=nT;
            }
            if (max==0) {
                consensus[pos]='-';
            }
        }
    }
}*/

void FindConsensus2 (string& consensus, vector<char>& alphabet, vector<string>& seqs) {
    consensus=seqs[0];
    if (seqs.size()>1) {
        for (int pos=0;pos<seqs[0].size();pos++) {
            std::unordered_map<char, int> counts;
            for (int seq=0;seq<seqs.size();seq++) {
                //Check if this is in the alphabet
                int check=CheckAlphabet(alphabet,seqs[seq][pos]);
                if (check==1) {
                    counts[seqs[seq][pos]]++;
                }
            }
            char most_common = 0;
            int max_count = 0;
            for (auto &[c, count] : counts) {
                if (count > max_count) {
                    max_count = count;
                    most_common = c;
                }
            }
            consensus[pos]=most_common;
        }
    }
}

int CheckAlphabet (vector<char>& alphabet, char a) {
    int match=0;
    auto it = find(alphabet.begin(), alphabet.end(), a);
    if (it != alphabet.end()) {
        match=1;
    }
    return match;
}


void FindSVariants (vector<sparseseq>& variants, string& consensus, vector<char>& alphabet, vector<string>& seqs) {
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<seqs[i].size();pos++) {
            if (seqs[i].compare(pos,1,consensus,pos,1)!=0) {
                int check=CheckAlphabet(alphabet,seqs[i][pos]);
                if (check==1) {
                    //cout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << seqs[i][pos] << "\n";
                    s.locus.push_back(pos);
                    s.allele.push_back(seqs[i][pos]);
                }
            }
        }
        variants.push_back(s);
    }
}

void FindPairwiseDistances (vector<sparseseq>& variants, vector<string>& seqs, vector<char>& alphabet, vector< vector<int> >& seqdists) {
    vector<int> zeros(seqs.size(),0);
    for (int i=0;i<seqs.size();i++) {
        seqdists.push_back(zeros);
    }
    for (int i=0;i<seqs.size();i++) {
        for (int j=i+1;j<seqs.size();j++) {
            int dist=0;
            //Find unique difference positions;
            vector<int> uniq;
            for (int k=0;k<variants[i].locus.size();k++) {
                uniq.push_back(variants[i].locus[k]);
            }
            for (int k=0;k<variants[j].locus.size();k++) {
                uniq.push_back(variants[j].locus[k]);
            }
            sort(uniq.begin(),uniq.end());
            uniq.erase(unique(uniq.begin(),uniq.end()),uniq.end());
            for (int k=0;k<uniq.size();k++) {
                int check1=CheckAlphabet(alphabet,seqs[i][uniq[k]]);
                int check2=CheckAlphabet(alphabet,seqs[j][uniq[k]]);
                if (check1==1&&check2==1) {
                    if (seqs[i][uniq[k]]!=seqs[j][uniq[k]]) {
                        dist++;
                    }
                }
            }
            seqdists[i][j]=dist;
            seqdists[j][i]=dist;
        }
    }
    for (int i=0;i<seqs.size();i++) {
        if (seqs[i].size()==0) {
            for (int j=0;j<seqs.size();j++){
                seqdists[i][j]=-1;
                seqdists[j][i]=-1;
            }
            seqdists[i][i]=0;
        }
    }
}


void GetSubsetsIJ (run_params p, const vector<string>& names, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets) {
    FindDistanceSubsetsIJ(p.dist_cut,seqdists,subsets);
    ofstream subs_file;
    subs_file.open("Subset_data.out");
    cout << "Subsets\n";
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            subs_file << names[subsets[i][j]] << " ";
        }
        subs_file << "\n";
    }
}


void FindDistanceSubsetsIJ(int cut, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets) {
    vector <int> range;
    for (int i=0;i<seqdists.size();i++) {
        range.push_back(i);
    }
    while (range.size()>0) {
        vector<int> sset;
        sset.clear();
        sset.push_back(range[0]);
        int add=1;
        while (add==1) {
            int ss=sset.size();
            add=0;
            for (int i=0;i<ss;i++) {
                for (int j=0;j<seqdists.size();j++) {
                    if (seqdists[sset[i]][j]<=cut) {
                        //cout << "Add " << sset[i] << " " << j << "\n";
                        sset.push_back(j);
                    }
                    if (seqdists[j][sset[i]]<=cut) {
                        //cout << "Add " << sset[i] << " " << j << "\n";
                        sset.push_back(j);
                    }
                }
            }
            sort(sset.begin(),sset.end());
            sset.erase(unique(sset.begin(),sset.end()),sset.end());
            //cout << sset.size() << "\n";
            if (sset.size()>ss) {
                add=1;
            }
        }
        //Add sset to subsets
        //cout << "Push back to subsets\n";
        subsets.push_back(sset);
        for (int i=0;i<sset.size();i++) {
            for (int j=0;j<range.size();j++) {
                if (range[j]==sset[i]) {
                    range.erase(range.begin()+j);
                    break;
                }
            }
        }
    }
}


/*void GetAliStats (const vector<string>& seqs,vector<site>& ali_stats) {
    for (int i=0;i<seqs.size();i++) {
        if (ali_stats.size()==0) {
            for (int j=0;j<seqs[i].length();j++) {
                site s;
                s.A=0;
                s.C=0;
                s.G=0;
                s.T=0;
                s.N=0;
                if (seqs[i].compare(j,1,"A")==0) {
                    s.A++;
                    s.N++;
                }
                if (seqs[i].compare(j,1,"C")==0) {
                    s.C++;
                    s.N++;
                }
                if (seqs[i].compare(j,1,"G")==0) {
                    s.G++;
                    s.N++;
                }
                if (seqs[i].compare(j,1,"T")==0) {
                    s.T++;
                    s.N++;
                }
                ali_stats.push_back(s);
            }
        } else {
            for (int j=0;j<seqs[i].length();j++) {
                if (seqs[i].compare(j,1,"A")==0) {
                    ali_stats[j].A++;
                    ali_stats[j].N++;
                }
                if (seqs[i].compare(j,1,"C")==0) {
                    ali_stats[j].C++;
                    ali_stats[j].N++;
                }
                if (seqs[i].compare(j,1,"G")==0) {
                    ali_stats[j].G++;
                    ali_stats[j].N++;
                }
                if (seqs[i].compare(j,1,"T")==0) {
                    ali_stats[j].T++;
                    ali_stats[j].N++;
                }
            }
        }
    }
}*/


void GetAliStats2 (run_params& p, const vector<string>& seqs, vector<char>& alphabet, vector<site2>& ali_stats) {
    if (p.verb==1) {
        cout << "Getting alignment statistics\n";
    }
    for (int i=0;i<seqs.size();i++) {
        if (ali_stats.size()==0) {
            for (int j=0;j<seqs[i].length();j++) {
                site2 s;
                s.N=0;
                for (int i=0;i<alphabet.size();i++) {
                    s.counts.push_back(0);
                }
                for (int k=0;k<alphabet.size();k++) {
                    if (seqs[i][j]==alphabet[k]) {
                        s.counts[k]++;
                        s.N++;
                    }
                }
                ali_stats.push_back(s);
            }
        } else {
            for (int j=0;j<seqs[i].length();j++) {
                for (int k=0;k<alphabet.size();k++) {
                    if (seqs[i][j]==alphabet[k]) {
                        ali_stats[j].counts[k]++;
                        ali_stats[j].N++;
                    }
                }
            }
        }
    }
}


void FindVariants (vector<site2>& ali_stats, vector<int>& var_positions) {
    for (int i=0;i<ali_stats.size();i++) {
        int count=0;
        for (int j=0;j<ali_stats[i].counts.size();j++) {
            if (ali_stats[i].counts[j]>0) {
                count++;
            }
        }
        if (count>1) {
            ali_stats[i].variant=1;
            var_positions.push_back(i);
        } else {
            ali_stats[i].variant=0;
        }
    }
}

/*void GetConsensus (vector<site>& ali_stats, vector<string>& consensus) {
    ofstream cons_file;
    cons_file.open("Alignment_consensus.fa");
    cons_file << ">Alignment_consensus\n";
    for (int i=0;i<ali_stats.size();i++) {
        string cons="N";
        if (ali_stats[i].A>0) {
            cons="A";
        }
        if (ali_stats[i].C>ali_stats[i].A) {
            cons="C";
        }
        if (ali_stats[i].G>ali_stats[i].A&&ali_stats[i].G>ali_stats[i].C) {
            cons="G";
        }
        if (ali_stats[i].T>ali_stats[i].G&&ali_stats[i].T>ali_stats[i].A&&ali_stats[i].T>ali_stats[i].C) {
            cons="T";
        }
        cons_file << cons;
        consensus.push_back(cons);
    }
    cons_file << "\n";
}*/

void GetConsensus2 (vector<site2>& ali_stats, vector<char>& alphabet, vector<string>& consensus) {
    ofstream cons_file;
    cons_file.open("Alignment_consensus.fa");
    cons_file << ">Alignment_consensus\n";
    for (int i=0;i<ali_stats.size();i++) {
        string cons="N";
        int max=0;
        for (int j=0;j<alphabet.size();j++) {
            if (ali_stats[i].counts[j]>max) {
                max=ali_stats[i].counts[j];
                cons=alphabet[j];
            }
        }
        cons_file << cons;
        consensus.push_back(cons);
    }
    cons_file << "\n";
}


void CalculateFrequencies (run_params& p, vector<site2>& ali_stats, vector<char>& alphabet, vector<string>& second) {
    int index=alphabet.size();
    if (p.verb==1) {
        cout << "Frequencies\n";
    }
    for (int i=0;i<ali_stats.size();i++) {
        if (ali_stats[i].variant==1) {
            vector<double> c;
            for (int j=0;j<ali_stats[i].counts.size();j++) {
                c.push_back(ali_stats[i].counts[j]);
            }
            sort(c.begin(),c.end()); //c[0] is the smallest
            for (int j=0;j<ali_stats[i].counts.size();j++) {
                if (c[index-2]==ali_stats[i].counts[j]) {
                    second[i]=alphabet[j];
                }
            }

            ali_stats[i].freq=c[index-2]/(c[index-2]+c[index-1]);
            if (p.verb==1) {
                cout << i << " ";
                for (int j=0;j<ali_stats[i].counts.size();j++) {
                    cout << ali_stats[i].counts[j] << " ";
                }
                cout << second[i] << " " << ali_stats[i].freq << "\n";
            }
        }
    }
}



void MakeInitialPairs (vector<int>& var_positions, vector<pr>& pairs) {
    for (int i=0;i<var_positions.size();i++) {
        for (int j=0;j<var_positions.size();j++) {
            pr p;
            p.i=var_positions[i];
            p.j=var_positions[j];
            p.c11=0;
            p.c10=0;
            p.c01=0;
            p.c00=0;
            pairs.push_back(p);
        }
    }
}

void ConstructPairs (run_params p, const vector<string>& second, vector<pr>& pairs) {
    ifstream ali_file;
    ali_file.open(p.ali_file);
    for (int i=0;i<100000;i++) {
        string name;
        if (!(ali_file >> name)) break;
        if (name.compare(0,1,">")!=0) {
            for (int j=0;j<pairs.size();j++) {
                int p1=0;
                int p2=0;
                if (name.compare(pairs[j].i,1,second[pairs[j].i])==0) {
                    p1=1;
                }
                if (name.compare(pairs[j].j,1,second[pairs[j].j])==0) {
                    p2=1;
                }
                if (p1==1&&p2==1) {
                    pairs[j].c11++;
                }
                if (p1==1&&p2==0) {
                    pairs[j].c10++;
                }
                if (p1==0&&p2==1) {
                    pairs[j].c01++;
                }
                if (p1==0&&p2==0) {
                    pairs[j].c00++;
                }
            }
        }
    }
}


void FindCorrelations (vector<site2>& ali_stats, vector<pr>& pairs) {
    cout << "Find correlations\n";
    for (int i=0;i<pairs.size();i++) {
        double top=pairs[i].c00*(-ali_stats[pairs[i].i].freq)*(-ali_stats[pairs[i].j].freq);
        top=top+pairs[i].c01*(-ali_stats[pairs[i].i].freq)*(1-ali_stats[pairs[i].j].freq);
        top=top+pairs[i].c10*(1-ali_stats[pairs[i].i].freq)*(-ali_stats[pairs[i].j].freq);
        top=top+pairs[i].c11*(1-ali_stats[pairs[i].i].freq)*(1-ali_stats[pairs[i].j].freq);
        double b1=0;
        b1=b1+(pairs[i].c00+pairs[i].c01)*pow(-ali_stats[pairs[i].i].freq,2);
        b1=b1+(pairs[i].c10+pairs[i].c11)*pow(1-ali_stats[pairs[i].i].freq,2);
        double b2=0;
        b2=b2+(pairs[i].c00+pairs[i].c10)*pow(-ali_stats[pairs[i].j].freq,2);
        b2=b2+(pairs[i].c01+pairs[i].c11)*pow(1-ali_stats[pairs[i].j].freq,2);
        double b=b1*b2;
        b=sqrt(b);
        pairs[i].correl=(top+0.)/(b+0.);
        cout << pairs[i].c11 << " " << pairs[i].c10 << " " << pairs[i].c01 << " " << pairs[i].c00 << " " << top << " " << b << " " << pairs[i].correl << "\n";
    }
}


void FindIdentical (const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident) {
    for (int i=0;i<var_positions.size();i++) {
        vector<int> id;
        id.push_back(i);
        for (int j=i+1;j<var_positions.size();j++) {
            int diff=0;
            for (int k=0;k<var_positions.size();k++) {
                if (correls[i][k]!=correls[j][k]) {
                    diff=1;
                    break;
                }
            }
            if (diff==0) {
                id.push_back(j);
            }
        }
        ident.push_back(id);
    }
}

void FindIDCorrel (double tol, const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident) {
    for (int i=0;i<var_positions.size();i++) {
        vector<int> id;
        id.push_back(i);
        for (int j=i+1;j<var_positions.size();j++) {
            if (correls[i][j]>1-tol) {
                id.push_back(j);
            }
        }
        ident.push_back(id);
    }
}


void FindRemovals (const vector< vector<int> >& ident, vector<int>& to_rem) {
    for (int i=0;i<ident.size();i++) {
        for (int j=1;j<ident[i].size();j++) {
            to_rem.push_back(ident[i][j]);
        }
    }
    sort(to_rem.begin(),to_rem.end());
    reverse(to_rem.begin(),to_rem.end());
    to_rem.erase(unique(to_rem.begin(),to_rem.end()),to_rem.end());

}

void DoRemovals (const vector<int>& to_rem, vector< vector<double> >& correls, vector< vector<int> >& ident, vector<double>& frequencies) {
    for (int i=0;i<to_rem.size();i++) {
        for (int j=0;j<correls.size();j++) {
            correls[j].erase(correls[j].begin()+to_rem[i]);
        }
    }
    for (int i=0;i<to_rem.size();i++) {
        correls.erase(correls.begin()+to_rem[i]);
        ident.erase(ident.begin()+to_rem[i]);
        frequencies.erase(frequencies.begin()+to_rem[i]);
    }
}


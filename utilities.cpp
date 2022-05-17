#include "aligntools.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void GetParameters (run_params& p, int argc, const char **argv) {
	string p_switch;
	int x=1;
	p.ali_file="Align.fa";
    p.get_positions=1;
    p.get_frequencies=1;
    p.get_correlations=1;
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
        } else if (p_switch.compare("--cutoff")==0) {
            x++;
            p.dist_cut=atoi(argv[x]);
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

void CheckBaseCase (vector<string>& seqs) {
    cout << "Check base case\n";
    for (int i=0;i<seqs.size();i++) {
        for (int j=0;j<seqs[i].size();j++) {
            if (seqs[i].compare(j,1,"a")==0) {
                seqs[i][j]='A';
            } else if (seqs[i].compare(j,1,"c")==0) {
                seqs[i][j]='C';
            } else if (seqs[i].compare(j,1,"g")==0) {
                seqs[i][j]='G';
            } else if (seqs[i].compare(j,1,"t")==0) {
                seqs[i][j]='T';
            } else if (seqs[i].compare(j,1,"n")==0) {
                seqs[i][j]='N';
            } else if (seqs[i].compare(j,1,"-")==0) {
                seqs[i][j]='N';
            }
        }
    }
}

void FindConsensus (string& consensus, vector<string>& seqs) {
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
}

void FindSVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs) {
    cout << "Find variants\n";
    for (int i=0;i<seqs.size();i++) {
        sparseseq s;
        for (int pos=0;pos<seqs[i].size();pos++) {
            if (seqs[i].compare(pos,1,consensus,pos,1)!=0) {
                if (seqs[i].compare(pos,1,"A")==0||seqs[i].compare(pos,1,"C")==0||seqs[i].compare(pos,1,"G")==0||seqs[i].compare(pos,1,"T")==0) {
                    //cout << "Found variant " << pdat[i].code_match << " " << pos << " " << consensus[pos] << " " << seqs[i][pos] << "\n";
                    s.locus.push_back(pos);
                    s.allele.push_back(seqs[i][pos]);
                }
            }
        }
        variants.push_back(s);
    }
}

void FindPairwiseDistances (vector<sparseseq>& variants, vector<string>& seqs, vector< vector<int> >& seqdists) {
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
                if (seqs[i][uniq[k]]=='A'||seqs[i][uniq[k]]=='C'||seqs[i][uniq[k]]=='G'||seqs[i][uniq[k]]=='T') {
                    if (seqs[j][uniq[k]]=='A'||seqs[j][uniq[k]]=='C'||seqs[j][uniq[k]]=='G'||seqs[j][uniq[k]]=='T') {
                        if (seqs[i][uniq[k]]!=seqs[j][uniq[k]]) {
                            dist++;
                        }
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
    cout << "Consistent groups\n";
    FindDistanceSubsetsIJ(p.dist_cut,seqdists,subsets);
    cout << "Subsets\n";
    for (int i=0;i<subsets.size();i++) {
        for (int j=0;j<subsets[i].size();j++) {
            cout << names[subsets[i][j]] << " ";
        }
        cout << "\n";
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




void ReadVariants (run_params& p, vector<site>& ali_stats) {
    ifstream ali_file;
    ali_file.open(p.ali_file);
    for (int i=0;i<1000000;i++) {
        string name;
        if (!(ali_file >> name)) break;
        if (name.compare(0,1,">")!=0) {
            if (ali_stats.size()==0) {
                for (int j=0;j<name.length();j++) {
                    site s;
                    s.A=0;
                    s.C=0;
                    s.G=0;
                    s.T=0;
                    if (name.compare(j,1,"a")==0||name.compare(j,1,"A")==0) {
                        s.A++;
                    }
                    if (name.compare(j,1,"c")==0||name.compare(j,1,"C")==0) {
                        s.C++;
                    }
                    if (name.compare(j,1,"g")==0||name.compare(j,1,"G")==0) {
                        s.G++;
                    }
                    if (name.compare(j,1,"t")==0||name.compare(j,1,"T")==0) {
                        s.T++;
                    }
                    ali_stats.push_back(s);
                }
            } else {
                for (int j=0;j<name.length();j++) {
                    if (name.compare(j,1,"a")==0||name.compare(j,1,"A")==0) {
                        ali_stats[j].A++;
                    }
                    if (name.compare(j,1,"c")==0||name.compare(j,1,"C")==0) {
                        ali_stats[j].C++;
                    }
                    if (name.compare(j,1,"g")==0||name.compare(j,1,"G")==0) {
                        ali_stats[j].G++;
                    }
                    if (name.compare(j,1,"t")==0||name.compare(j,1,"T")==0) {
                        ali_stats[j].T++;
                    }
                }
            }
//        } else {
//	  	cout << name << "\n";
	}
    }
    ali_file.close();
}


void FindVariants (vector<site>& ali_stats, vector<int>& var_positions) {
    for (int i=0;i<ali_stats.size();i++) {
        int count=0;
        if (ali_stats[i].A>0) {
            count++;
        }
        if (ali_stats[i].C>0) {
            count++;
        }
        if (ali_stats[i].G>0) {
            count++;
        }
        if (ali_stats[i].T>0) {
            count++;
        }
        if (count>1) {
            ali_stats[i].variant=1;
            var_positions.push_back(i);
        } else {
            ali_stats[i].variant=0;
        }
    }
}

void GetConsensus (vector<site>& ali_stats, vector<string>& consensus) {
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
}

void CalculateFrequencies (vector<site>& ali_stats, vector<string>& second) {
    for (int i=0;i<ali_stats.size();i++) {
        if (ali_stats[i].variant==1) {
            vector<double> c;
            c.push_back(ali_stats[i].A);
            c.push_back(ali_stats[i].C);
            c.push_back(ali_stats[i].G);
            c.push_back(ali_stats[i].T);
            sort(c.begin(),c.end());
            if (c[2]==ali_stats[i].A) {
                second[i]="A";
            }
            if (c[2]==ali_stats[i].C) {
                second[i]="C";
            }
            if (c[2]==ali_stats[i].G) {
                second[i]="G";
            }
            if (c[2]==ali_stats[i].T) {
                second[i]="T";
            }

            ali_stats[i].freq=c[2]/(c[2]+c[3]);
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


void FindCorrelations (vector<site>& ali_stats, vector<pr>& pairs) {
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
    }
}

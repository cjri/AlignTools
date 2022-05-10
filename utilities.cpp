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
        } else {
			cout << "Incorrect usage\n ";
			exit(1);
		}
		p_switch.clear();
		x++;
	}
    if (p.get_frequencies==0) {
        p.get_correlations=0;
    }
}

void ReadVariants (run_params& p, vector<site>& ali_stats) {
    ifstream ali_file;
    ali_file.open(p.ali_file);
    for (int i=0;i<100000;i++) {
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

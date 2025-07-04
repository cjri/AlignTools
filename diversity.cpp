#include "aligntools.h"
#include "distance_matrix.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void GetPiDiversity (run_params& p, int seq_length, vector<string>& seqs, vector<string>& names) {
        vector< vector<int> > seqdists;
        MakeDistanceMatrix (p,0,seqs,names,seqdists);
        double tot=0;
        double c=0;
        for (int i=0;i<seqdists.size();i++) {
            for (int j=0;j<i;j++) {
                tot=tot+seqdists[i][j];
                c++;
            }
        }
        tot=tot/(seq_length+0.);
        tot=tot/c;
        cout << "Pi = " << tot << "\n";
}



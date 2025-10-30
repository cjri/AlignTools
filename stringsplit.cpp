#include "aligntools.h"
#include "distance_matrix.h"
#include "io.h"
#include "stringsplit.h"
#include "timesplit.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void SSplit (run_params& p, vector<string>& consensus, vector<int>& var_positions, vector<site2>& ali_stats, vector<string>& names, vector<string>& seqs) {
    vector<string> strings;
    vector<string> strings_uniq;
    ReadStrings(strings,strings_uniq);

    //For each time, generate a sub-alignment of the sequences
    //Time-dependent sequences
    vector< vector<string> > seqs_t;
    vector< vector<string> > names_t;
    SplitSeqsNames (strings,strings_uniq,seqs,names,seqs_t,names_t);
    
    //Output sequences
    for (int i=0;i<strings_uniq.size();i++) {
        string name=modifyFilenameS(p.ali_file,strings_uniq[i]);
        ofstream time_file;
        time_file.open(name.c_str());
        for (int j=0;j<seqs_t[i].size();j++) {
            time_file << ">" << names_t[i][j] << "\n";
            time_file << seqs_t[i][j] << "\n";
        }
        time_file.close();
    }
}

void SplitSeqsNames (const vector<string>& strings, const vector<string>& strings_uniq, const vector<string>& seqs, const vector<string>& names, vector< vector<string> >& seqs_t, vector< vector<string> >& names_t) {
    for (int i=0;i<strings.size();i++) {
        vector<string> s;
        vector<string> n;
        for (int j=0;j<seqs.size();j++) {
            if (strings[j]==strings_uniq[i]) {
                s.push_back(seqs[j]);
                n.push_back(names[j]);
            }
        }
        seqs_t.push_back(s);
        names_t.push_back(n);
    }
}

string modifyFilenameS(string ali_file, string s) {
    // Find the last dot in the filename
    size_t dotPos = ali_file.find_last_of('.');
    if (dotPos == string::npos) {
        return ali_file + s;
    }
    string start = ali_file.substr(0, dotPos);
    string end  = ali_file.substr(dotPos);
    return start +  "_" + s + end;
}


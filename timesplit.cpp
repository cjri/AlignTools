#include "aligntools.h"
#include "distance_matrix.h"
#include "io.h"
#include "timesplit.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

void GetTDNucleotideCounts (run_params& p, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& seqs) {
    vector<int> times;
    vector<int> times_uniq;
    ReadTimes(times,times_uniq);

    //For each time, generate a sub-alignment of the sequences
    //Time-dependent sequences
    vector< vector<string> > seqs_t;
    SplitSeqs (times,times_uniq,seqs,seqs_t);
    
    //Alignment statistics
    vector< vector<site> > ali_stats_t;
    SplitAliStats (times_uniq,seqs_t,ali_stats_t);
    
    //Time-dependent variant frequencies
    GetTimeVariants (times_uniq,ali_stats_t);
    
    vector<string> second=consensus;
    CalculateFrequencies (ali_stats,second);

    for (int i=0;i<times_uniq.size();i++) {
        vector<string> temp=consensus;
        for (int j=0;j<ali_stats_t[i].size();j++) {
            ali_stats_t[i][j].freq=0;
        }
        CalculateFrequencies (ali_stats_t[i],temp);
    }
    
    //Go through the frequencies of each variant with time
    OutputNucleotideCountsTime (p,consensus,second,times_uniq,var_positions,ali_stats_t);

}

void SplitSeqs (const vector<int>& times, const vector<int>& times_uniq, const vector<string>& seqs, vector< vector<string> >& seqs_t) {
    for (int i=0;i<times_uniq.size();i++) {
        vector<string> s;
        for (int j=0;j<seqs.size();j++) {
            if (times[j]==times_uniq[i]) {
                s.push_back(seqs[j]);
            }
        }
        seqs_t.push_back(s);
    }
}

void SplitAliStats (const vector<int>& times_uniq, const vector< vector<string> >& seqs_t, vector< vector<site> >& ali_stats_t) {
    for (int i=0;i<times_uniq.size();i++) {
        vector<site> a;
        GetAliStats (seqs_t[i],a);
        ali_stats_t.push_back(a);
    }
}

void GetTimeVariants (const vector<int>& times_uniq, vector< vector<site> >& ali_stats_t) {
    for (int i=0;i<times_uniq.size();i++) {
        vector<int> v;
        FindVariants (ali_stats_t[i],v);
    }
}

void TDSplit (run_params& p, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& names, vector<string>& seqs) {
    vector<int> times;
    vector<int> times_uniq;
    ReadTimes(times,times_uniq);

    //For each time, generate a sub-alignment of the sequences
    //Time-dependent sequences
    vector< vector<string> > seqs_t;
    vector< vector<string> > names_t;
    SplitSeqsNames (times,times_uniq,seqs,names,seqs_t,names_t);
    
    //Output sequences
    for (int i=0;i<times_uniq.size();i++) {
        string name=modifyFilename(p.ali_file,times_uniq[i]);
        ofstream time_file;
        time_file.open(name.c_str());
        for (int j=0;j<seqs_t[i].size();j++) {
            time_file << names_t[i][j] << "\n";
            time_file << seqs_t[i][j] << "\n";
        }
        time_file.close();
        
        
    }
}

void SplitSeqsNames (const vector<int>& times, const vector<int>& times_uniq, const vector<string>& seqs, const vector<string>& names, vector< vector<string> >& seqs_t, vector< vector<string> >& names_t) {
    for (int i=0;i<times_uniq.size();i++) {
        vector<string> s;
        vector<string> n;
        for (int j=0;j<seqs.size();j++) {
            if (times[j]==times_uniq[i]) {
                s.push_back(seqs[j]);
                n.push_back(names[j]);
            }
        }
        seqs_t.push_back(s);
        names_t.push_back(n);
    }
}


string modifyFilename(string ali_file, int i) {
    // Find the last dot in the filename
    size_t dotPos = ali_file.find_last_of('.');
    if (dotPos == string::npos) {
        return ali_file + to_string(i);
    }
    string start = ali_file.substr(0, dotPos);
    string end  = ali_file.substr(dotPos);
    return start +  "_" + to_string(i) + end;
}

void GetTDNucleotideCounts (run_params& p, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& seqs);
void SplitSeqs (const vector<int>& times, const vector<int>& times_uniq, const vector<string>& seqs, vector< vector<string> >& seqs_t);
void SplitAliStats (const vector<int>& times_uniq, const vector< vector<string> >& seqs_t, vector< vector<site> >& ali_stats_t);
void GetTimeVariants (const vector<int>& times_uniq, vector< vector<site> >& ali_stats_t);
void TDSplit (run_params& p, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& names, vector<string>& seqs);
void SplitSeqsNames (const vector<int>& times, const vector<int>& times_uniq, const vector<string>& seqs, const vector<string>& names, vector< vector<string> >& seqs_t, vector< vector<string> >& names_t);
string modifyFilename(string ali_file, int i);


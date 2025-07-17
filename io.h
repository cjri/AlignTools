
void ReadFastaAli (run_params p, vector<string>& seqs, vector<string>& names);
void ReadTimes (vector<int>& times, vector<int>& times_uniq);
void PrintVariantPositions (vector<int>& var_positions, vector<string>& consensus, vector<string>& second);
void PrintFrequencies (vector<double>& frequencies);
void PrintCorrelations (vector< vector<double> >& correls);
void OutputBitstrings (run_params& p, vector< vector<int> >& bitstrings, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<delet> deletions, vector< vector<int> >& denovo);
void OutputBitstringsBinary (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void GetDenovoSites (vector< vector<int> >& denovo, vector<int>& all_denovo);
void OutputSites (const vector<int>& var_positions, vector<int>& all_denovo, vector< vector<int> >& bitstrings);

void OutputBitstringsSparse (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void OutputBitstringsFasta (const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void EditDenovo (vector< vector<int> >& denovo);
void OutputNucleotideCountsTime (run_params& p, vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site> >& ali_stats_t);
void OutputVariantFrequencies (vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site> >& ali_stats_t);

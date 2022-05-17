
void GetParameters (run_params& p, int argc, const char **argv);
void ReadFastaAli (run_params p, vector<string>& seqs, vector<string>& names);
void CheckBaseCase (vector<string>& seqs);
void FindConsensus (string& consensus, vector<string>& seqs);
void FindSVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindPairwiseDistances (vector<sparseseq>& variants, vector<string>& seqs, vector< vector<int> >& seqdists);
void GetSubsetsIJ (run_params p, const vector<string>& names, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);
void FindDistanceSubsetsIJ(int cut, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);


void ReadVariants (run_params& p, vector<site>& ali_stats);
void FindVariants (vector<site>& ali_stats, vector<int>& var_positions);
void GetConsensus (vector<site>& ali_stats, vector<string>& consensus);
void CalculateFrequencies (vector<site>& ali_stats, vector<string>& second);
void MakeInitialPairs (vector<int>& var_positions, vector<pr>& pairs);
void ConstructPairs (run_params p, const vector<string>& second, vector<pr>& pairs);
void FindCorrelations (vector<site>& ali_stats, vector<pr>& pairs);


void GetParameters (run_params& p, int argc, const char **argv);
void CheckBaseCase (vector<string>& seqs);
void FindConsensus (string& consensus, vector<string>& seqs);
void FindSVariants (vector<sparseseq>& variants, string& consensus, vector<string>& seqs);
void FindPairwiseDistances (vector<sparseseq>& variants, vector<string>& seqs, vector< vector<int> >& seqdists);
void GetSubsetsIJ (run_params p, const vector<string>& names, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);
void FindDistanceSubsetsIJ(int cut, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);



void GetAliStats (const vector<string>& seqs,vector<site>& ali_stats);


void FindVariants (vector<site>& ali_stats, vector<int>& var_positions);
void GetConsensus (vector<site>& ali_stats, vector<string>& consensus);
void CalculateFrequencies (vector<site>& ali_stats, vector<string>& second);

void MakeInitialPairs (vector<int>& var_positions, vector<pr>& pairs);
void ConstructPairs (run_params p, const vector<string>& second, vector<pr>& pairs);
void FindCorrelations (vector<site>& ali_stats, vector<pr>& pairs);

void FindIdentical (const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident);

void FindIDCorrel (double tol, const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident);


void FindRemovals (const vector< vector<int> >& ident, vector<int>& to_rem);
void DoRemovals (const vector<int>& to_rem, vector< vector<double> >& correls, vector< vector<int> >& ident, vector<double>& frequencies);



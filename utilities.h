void GetParameters (run_params& p, int argc, const char **argv);
void CheckBaseCase (vector<string>& seqs);
void GetAlignmentType (run_params& p, vector<char>& alphabet, vector<string>& seqs);

void FindConsensus (string& consensus, vector<string>& seqs);
void FindConsensus2 (string& consensus, vector<char>& alphabet, vector<string>& seqs);
int CheckAlphabet (vector<char>& alphabet, char a);

void FindSVariants (vector<sparseseq>& variants, string& consensus, vector<char>& alphabet, vector<string>& seqs);
void FindPairwiseDistances (vector<sparseseq>& variants, vector<string>& seqs, vector<char>& alphabet, vector< vector<int> >& seqdists);
void GetSubsetsIJ (run_params p, const vector<string>& names, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);
void FindDistanceSubsetsIJ(int cut, const vector< vector<int> >& seqdists, vector< vector<int> >& subsets);



//void GetAliStats (const vector<string>& seqs,vector<site>& ali_stats);
void GetAliStats2 (run_params& p, const vector<string>& seqs, vector<char>& alphabet, vector<site2>& ali_stats);


void FindVariants (vector<site2>& ali_stats, vector<int>& var_positions);
//void GetConsensus (vector<site>& ali_stats, vector<string>& consensus);
void GetConsensus2 (vector<site2>& ali_stats, vector<char>& alphabet, vector<string>& consensus);

void CalculateFrequencies (run_params& p, vector<site2>& ali_stats, vector<char>& alphabet, vector<string>& second);

void MakeInitialPairs (vector<int>& var_positions, vector<pr>& pairs);
void ConstructPairs (run_params p, const vector<string>& second, vector<pr>& pairs);
void FindCorrelations (vector<site2>& ali_stats, vector<pr>& pairs);

void FindIdentical (const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident);

void FindIDCorrel (double tol, const vector<int>& var_positions, const vector< vector<double> >& correls, vector< vector<int> >& ident);


void FindRemovals (const vector< vector<int> >& ident, vector<int>& to_rem);
void DoRemovals (const vector<int>& to_rem, vector< vector<double> >& correls, vector< vector<int> >& ident, vector<double>& frequencies);



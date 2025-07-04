
void GenerateRandomSequences (run_params& p, int seq_length, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& seqs);

void ExtractFrequencies (run_params& p, vector<int>& var_positions, vector<site>& ali_stats, vector<double>& frequencies);
void ExtractCorrelations (run_params& p, vector<int>& var_positions, vector<string>& second, vector<site>& ali_stats, vector<pr>& pairs, vector< vector<double> >& correls);
void GenerateBitstrings (run_params& p, int n, int seq_length, const vector<string>& consensus, const vector<string>& second, const vector<int>& var_positions, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, vector<double>& var_freqs_cut, vector<double>& var_freqs_rem, vector< vector<double> >& correls_cut, vector<string>& seqs);
void UncorrelatedCorrection (run_params& p, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, const vector<double>& var_freqs_cut, const vector<double>& var_freqs_rem, vector< vector<int> >& bitstrings_cut, vector< vector<int> >& bitstrings);
void GetDeletions (vector<string>& seqs, vector<delet>& deletions);
vector< vector<int> > RandomString (int dim, int seq_length, int n, double& tolerance, const Eigen::VectorXd& mu, const Eigen::MatrixXd& omega, const Eigen::VectorXi& positions);



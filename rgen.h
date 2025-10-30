
void GenerateRandomSequences (run_params& p, int seq_length, vector<string>& consensus, vector<char>& alphabet, vector<int>& var_positions, vector<site2>& ali_stats, vector<string>& seqs, gsl_rng *rgen);

void ExtractFrequencies (run_params& p, vector<int>& var_positions, vector<site2>& ali_stats, vector<double>& frequencies);
void ExtractCorrelationsDels (run_params& p, vector<int>& conf_variants, vector<delet> deletions, vector<string>& second, vector<site2>& ali_stats, vector<pr>& pairs, vector<string>& seqs, vector< vector<double> >& correls);
void MakeInitialPairsDels (vector<int>& conf_variants, vector<delet>& deletions, vector<pr>& pairs);
void ConstructPairsDels (const vector<string>& second, vector<string> seqs, vector<delet>& deletions, vector<pr>& pairs);
void FindCorrelationsDels (run_params& p, vector<double>& qi, vector<double>& qj, vector<pr>& pairs);




void ExtractCorrelations (run_params& p, vector<int>& var_positions, vector<string>& second, vector<site2>& ali_stats, vector<pr>& pairs, vector< vector<double> >& correls);

void GetDenovo (vector<string>& consensus, vector<int>& delsites, vector<int>& conf_variants, vector< vector<int> >& bitstrings, vector<string>& seqs, vector< vector<int> >& denovo, gsl_rng *rgen);



void GenerateBitstrings (run_params& p, int n, int seq_length, const vector<string>& consensus, const vector<string>& second, const vector<int>& var_positions, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, vector<double>& var_freqs_cut, vector<double>& var_freqs_rem, vector< vector<double> >& correls_cut, vector<string>& seqs, gsl_rng *rgen);
void UncorrelatedCorrection (run_params& p, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, const vector<double>& var_freqs_cut, const vector<double>& var_freqs_rem, vector< vector<int> >& bitstrings_cut, vector< vector<int> >& bitstrings);

void GetDeletions (vector<string>& seqs, vector<delet>& deletions);
void FindDeletedSites50 (const vector<delet>& deletions, vector<int>& delsites);
void GetConfVariants (const vector<int>& var_positions, const vector<double>& frequencies, vector<int> delsites, vector<int>& conf_variants, vector<double>& conf_frequencies);
void GetInvariantSites (const vector<string>& consensus, vector<int>& conf_variants, vector<int>& delsites, vector<int>& inv_sites);

vector< vector<int> > RandomString (int dim, int seq_length, int n, double& tolerance, const Eigen::VectorXd& mu, const Eigen::MatrixXd& omega, const Eigen::VectorXi& positions);
double KingmanRatio (int n);



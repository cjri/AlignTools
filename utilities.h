
void GetParameters (run_params& p, int argc, const char **argv);
void ReadVariants (run_params& p, vector<site>& ali_stats);
void FindVariants (vector<site>& ali_stats, vector<int>& var_positions);
void GetConsensus (vector<site>& ali_stats, vector<string>& consensus);
void CalculateFrequencies (vector<site>& ali_stats, vector<string>& second);
void MakeInitialPairs (vector<int>& var_positions, vector<pr>& pairs);
void ConstructPairs (run_params p, const vector<string>& second, vector<pr>& pairs);
void FindCorrelations (vector<site>& ali_stats, vector<pr>& pairs);


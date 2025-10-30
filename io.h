
void ReadFastaAli (run_params& p, vector<string>& seqs, vector<string>& names);
void ReadTimes (vector<int>& times, vector<int>& times_uniq);
void ReadStrings (vector<string>& strings, vector<string>& strings_unique);

void PrintVariantPositions (vector<int>& var_positions, vector<string>& consensus, vector<string>& second);
void PrintFrequencies (vector<double>& frequencies);
void PrintCorrelations (vector< vector<double> >& correls);
void OutputBitstrings (run_params& p, vector< vector<int> >& bitstrings, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<char>& alphabet, vector<int>& conf_variants, vector<delet>& deletions, vector< vector<int> >& denovo);
void GetActualDeletions (run_params& p, vector<int>& conf_variants, vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& actual_del);
void MakeBitstringsNoDel (vector<int>& conf_variants, vector< vector<int> >& bitstrings);


void OutputBitstringsFasta (run_params& p, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<char>& alphabet, const vector<delet>& deletions, vector< vector<int> >& actual_del, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void PrintNextCharacter (int& indexs, const vector<string>& consensus, vector<char>& alphabet);
void OutputBitstringsSparse (run_params& p, const vector<int>& var_positions, const vector<string>& consensus, const vector<string>& second, vector<char>& alphabet, const vector<delet>& deletions, vector< vector<int> >& actual_del, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);




//void OutputBitstringsBinary (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void GetDenovoSites (vector< vector<int> >& denovo, vector<int>& all_denovo);
//void OutputSites (const vector<int>& var_positions, vector<int>& all_denovo, vector< vector<int> >& bitstrings);

//void OutputBitstringsSparse (const vector<int>& var_positions, const vector<delet>& deletions, vector< vector<int> >& bitstrings, vector< vector<int> >& denovo);
void EditDenovo (vector< vector<int> >& denovo);
void OutputNucleotideCountsTime (run_params& p, vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site2> >& ali_stats_t);
void OutputVariantFrequencies (vector<string>& consensus, vector<string>& second, const vector<int> times_uniq, const vector<int>& var_positions, const vector< vector<site2> >& ali_stats_t);
void OutputAlignmentFiltered (vector<string>& names, vector<string>& seqs, vector<int>& keep);
void OutputAlignmentSFiltered (vector<string>& names, vector<string>& seqs, vector<int>& sample);


void PrintDeletions (vector<delet>& deletions);
void PrintDelSites (vector<int>& delsites);
void PrintConfVar (vector<int>& conf_variants, vector<double>& conf_frequencies);
void PrintBitstrings (vector<int>& conf_variants, vector< vector<int> >& bitstrings);




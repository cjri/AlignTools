
#include <iostream>
#include <iomanip>
#include <fstream>
#include <cstring>
#include <vector>
#include <string>
#include <algorithm>
#include <map>
#include <cmath>

using namespace std;


struct run_params {
    int get_positions;
    int get_frequencies;
    int get_correlations;
    string ali_file;
    int dismat;
    int dist_cut;
};

struct sparseseq {
    vector<int> locus;
    vector<char> allele;
};

struct site {
    int A;
    int C;
    int G;
    int T;
    int variant; //Flag for variation at this site
    double freq; //Frequency of minor variant
};

struct pr {
    int i;
    int j;
    int c00;
    int c01;
    int c10;
    int c11;
    double correl;
};

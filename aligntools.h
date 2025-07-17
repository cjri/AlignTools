
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

#include <gsl/gsl_randist.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_eigen.h>
#include <gsl/gsl_cdf.h>
#include <gsl/gsl_rng.h>

#include <Eigen/Dense>


struct run_params {
    string method;
    int seed; //For random number generator
    int get_positions;
    int get_frequencies;
    int get_correlations;
    string ali_file;
    int dismat;
    int dist_cut;
    int n_generations;
    int verb;
    string output;
    double cutoff;
    double qq_cut; //Frequency cutoff to print variants
    int n_cut; //Depth cutoff to print variants
    int n_reps; //Number of times variant observed
    int denovo; //Add de novo variants to randomly-generated sequences
};

struct delet {
    int start;
    int length;
    double freq;
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
    int N;
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

#include <Eigen/Dense>
#include <random>
#include "aligntools.h"
#include "io.h"
#include "rgen.h"
#include "utilities.h"
#include <iostream>
#include <string>
#include <cstring>

using namespace Eigen;

void GenerateRandomSequences (run_params& p, int seq_length, vector<string>& consensus, vector<int>& var_positions, vector<site>& ali_stats, vector<string>& seqs,gsl_rng *rgen) {
    vector<string> second=consensus;
    CalculateFrequencies (ali_stats,second);
    
    if (p.verb==1) {
        PrintVariantPositions (var_positions,consensus,second);
    }

    //Get frequencies - binary
    vector<double> frequencies;
    ExtractFrequencies (p,var_positions,ali_stats,frequencies);
    
    vector<int> var_positions_cut;
    vector<int> var_positions_rem;
    vector<double> var_freqs_cut;
    vector<double> var_freqs_rem;
    if (p.cutoff>0) {
        for (int i=0;i<frequencies.size();i++) {
            if (frequencies[i]>p.cutoff) {
                var_positions_cut.push_back(var_positions[i]);
                var_freqs_cut.push_back(frequencies[i]);
            } else {
                var_positions_rem.push_back(var_positions[i]);
                var_freqs_rem.push_back(frequencies[i]);
            }
        }
    } else {
        var_positions_cut=var_positions;
        var_freqs_cut=frequencies;
    }
    
    vector<pr> pairs;
    vector< vector<double> > correls_cut;
    ExtractCorrelations (p,var_positions,second,ali_stats,pairs,correls_cut);
    
    GenerateBitstrings (p,p.n_generations,seq_length,consensus,second,var_positions,var_positions_cut,var_positions_rem,var_freqs_cut,var_freqs_rem,correls_cut,seqs,rgen);

}

void ExtractFrequencies (run_params& p, vector<int>& var_positions, vector<site>& ali_stats, vector<double>& frequencies) {
    if (p.get_frequencies==1) {
        for (int i=0;i<var_positions.size();i++) {
            frequencies.push_back(ali_stats[var_positions[i]].freq);
        }
        if (p.verb==1) {
            PrintFrequencies (frequencies);
        }
    }
}

void ExtractCorrelations (run_params& p, vector<int>& var_positions, vector<string>& second, vector<site>& ali_stats, vector<pr>& pairs, vector< vector<double> >& correls) {
    if (p.get_correlations==1) {
        MakeInitialPairs (var_positions,pairs);
        ConstructPairs (p,second,pairs);
        //Find correlations
        FindCorrelations (ali_stats,pairs);
        int index=0;
        for (int i=0;i<var_positions.size();i++) {
            vector<double> c;
            for (int j=0;j<var_positions.size();j++) {
                c.push_back(pairs[index].correl);
                index++;
            }
            correls.push_back(c);
        }
        if (p.verb==1) {
            PrintCorrelations (correls);
        }
    }
}

void GenerateBitstrings (run_params& p, int n, int seq_length, const vector<string>& consensus, const vector<string>& second, const vector<int>& var_positions, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, vector<double>& var_freqs_cut, vector<double>& var_freqs_rem, vector< vector<double> >& correls_cut, vector<string>& seqs, gsl_rng *rgen) {
    int dim=var_freqs_cut.size();
    
    //Eigen::VectorXd mu = Eigen::Map<Eigen::VectorXd>(frequencies.data(), frequencies.size())
    
    VectorXd mu(dim);
    for (int i=0;i<var_freqs_cut.size();i++) {
        mu(i) = var_freqs_cut[i];
    }
    
    MatrixXd omega(dim, dim);
    for (int i=0;i<dim; i++) {
        for (int j=0; j<dim; j++) {
            omega(i, j) = correls_cut[i][j];
        }
    }
    
    VectorXi positions(dim);
    for (int i=0;i<var_positions_cut.size();i++) {
        positions(i) = var_positions_cut[i];
    }
        
    //Find all deletions: Approach to handling indels
    vector<delet> deletions;
    GetDeletions (seqs,deletions);
    
    vector< vector<int> > bitstrings;
    vector< vector<int> > denovo;

    if (p.denovo==1) {
        //Find deletion sites of at least 50%
        vector<int> delsites;
        FindDeletedSites50(deletions,delsites);
        
        //Find variant positions that don't occur at deleted sites
        vector<int> conf_variants;
        GetConfVariants (var_positions,delsites,conf_variants);
        
        //Find invariant positions
        vector<int> inv_sites;
        GetInvariantSites (consensus,conf_variants,delsites,inv_sites);
        
        //invprop is the proportion of non-deleted sites with no variation
        double invprop=1-((conf_variants.size())/(consensus.size() - delsites.size()+0.));
        double invtot=consensus.size()-(conf_variants.size()+delsites.size());
        //cout << "Invariant proporion " << invprop << "\n";
        
        //Have proportion of sites with variants.
        
        
        
        double rate_denovo=0;
        if (invprop>0) {
            //Estimate rate of variant gain per site
            //This assumes rate is constant at all sites.
            //Calculated from the number of sites with no variants via a Poisson formula.
            
            double rate = -log(invprop);
            //cout << "Rate " << rate << "\n";
            
            //Estimate proportion of ratetime in new versus old tree
            double rat=KingmanRatio(seqs.size());
            
            //Assume evolution occurs at the given rate for the new branch
            rate_denovo=rate*rat*invtot;
            cout << "Here " << rate_denovo << "\n";
                        
        }
        
        double tolerance = 1e-6;
        //Call string generation
        vector< vector<int> > bitstrings_cut=RandomString (dim,seq_length,n,tolerance,mu,omega,positions);
        
        
        //Here add in the other frequencies
        UncorrelatedCorrection (p,var_positions_cut,var_positions_rem,var_freqs_cut,var_freqs_rem,bitstrings_cut,bitstrings);

        
        //Find positions of de novo mutations for each sequence
        for (int i=0;i<bitstrings.size();i++) {
            vector<int> dn;
            //double r=rate_denovo;
            //Numbers of de novo sites in each genome
            int ndn=gsl_ran_poisson(rgen,rate_denovo);
            vector<int> inv_sites_found;
            for (int j=0;j<inv_sites.size();j++) {
                inv_sites_found.push_back(0);
            }
            for (int j=0;j<ndn;j++) {
                int site=floor(gsl_rng_uniform(rgen)*inv_sites.size()+0.5)-1;
                if (inv_sites_found[site]==0) {
                    site=inv_sites[site];
                    dn.push_back(site);
                    inv_sites_found[site]=1;
                }
            }
            denovo.push_back(dn);
        }
        
        for (int i=0;i<denovo.size();i++) {
            cout << "String " << i << " " << rate_denovo << " ";
            for (int j=0;j<denovo[i].size();j++) {
                cout << denovo[i][j] << " ";
            }
            cout << "\n";
        }

        
    } else {
        double tolerance = 1e-6;
        //Call string generation
        vector< vector<int> > bitstrings_cut=RandomString (dim,seq_length,n,tolerance,mu,omega,positions);
        
        
        //Here add in the other frequencies
        UncorrelatedCorrection (p,var_positions_cut,var_positions_rem,var_freqs_cut,var_freqs_rem,bitstrings_cut,bitstrings);
    
    }
    
/*    cout << "Deletions\n";
    for (int k=0;k<deletions.size();k++) {
        cout << deletions[k].start << " " << deletions[k].length << " " << deletions[k].freq << "\n";
    }*/
    
    
    
    
    OutputBitstrings (p,bitstrings,var_positions,consensus,second,deletions,denovo);
        
}

void GetConfVariants (const vector<int>& var_positions, vector<int> delsites, vector<int>& conf_variants) {
    for (int i=0;i<var_positions.size();i++) {
        int deleted=0;
        for (int j=0;j<delsites.size();j++) {
            if (var_positions[i]==delsites[j]) {
                deleted=1;
                break;
            }
        }
        if (deleted==0) {
            conf_variants.push_back(var_positions[i]);
        }
    }
}

void GetInvariantSites (const vector<string>& consensus, vector<int>& conf_variants, vector<int>& delsites, vector<int>& inv_sites) {
    for (int i=0;i<consensus.size();i++) {
        int found=0;
        for (int j=0;j<conf_variants.size();j++) {
            if (conf_variants[j]==i) {
                found=1;
                break;
            }
        }
        if (found==0) {
            for (int j=0;j<delsites.size();j++) {
                if (delsites[j]==i) {
                    found=1;
                    break;
                }
            }
        }
        if (found==0) {
            inv_sites.push_back(i);
        }
    }
}

void GetDeletions (vector<string>& seqs, vector<delet>& deletions) {
    for (int i=0;i<seqs.size();i++) {
        int found=0;
        delet d;
        d.start=-1;
        d.length=0;
        d.freq=0;
        for (int j=0;j<seqs[i].size();j++) {
            if (seqs[i].compare(j,1,"N")!=0) {
                if (found==1) {
                    int seen=0;
                    for (int k=0;k<deletions.size();k++) {
                        if (deletions[k].start==d.start&&deletions[k].length==d.length) {
                            deletions[k].freq++;
                            seen=1;
                        }
                    }
                    if (seen==0) {
                        deletions.push_back(d);
                    }
                    found=0;
                    d.start=-1;
                    d.length=0;
                    d.freq=0;
                }
            }
            if (seqs[i].compare(j,1,"N")==0&&found==1) {
                d.length++;
            }
            if (seqs[i].compare(j,1,"N")==0&&found==0) {
                found=1;
                d.start=j;
                d.length=1;
                d.freq=1;
            }
        }
        if (found==1) {
            int seen=0;
            for (int k=0;k<deletions.size();k++) {
                if (deletions[k].start==d.start&&deletions[k].length==d.length) {
                    deletions[k].freq++;
                    seen=1;
                }
            }
            if (seen==0) {
                deletions.push_back(d);
            }
            
        }
    }
    for (int i=0;i<deletions.size();i++) {
        deletions[i].freq=deletions[i].freq/(seqs.size()+0.);
    }
}

void FindDeletedSites50 (const vector<delet>& deletions, vector<int>& delsites) {
    for (int i=0;i<deletions.size();i++) {
        if (deletions[i].freq>0.5) {
            for (int j=0;j<deletions[i].length;j++) {
                delsites.push_back(deletions[i].start+j);
            }
        }
    }
    sort(delsites.begin(),delsites.end());
    delsites.erase(unique(delsites.begin(),delsites.end()),delsites.end());
}


void UncorrelatedCorrection (run_params& p, const vector<int>& var_positions_cut, const vector<int>& var_positions_rem, const vector<double>& var_freqs_cut, const vector<double>& var_freqs_rem, vector< vector<int> >& bitstrings_cut, vector< vector<int> >& bitstrings) {
    if (p.cutoff>0) {
        for (int i=0;i<bitstrings_cut.size();i++) {
            vector<int> b;
            int index1=0;
            int index2=0;
            while (index1<var_positions_cut.size()&&index2<var_positions_rem.size()) {
                if (var_positions_cut[index1]<var_positions_rem[index2]) {
                    b.push_back(bitstrings_cut[i][index1]);
                    index1++;
                } else {
                    random_device rd;
                    mt19937 gen(rd());
                    uniform_real_distribution<> dis(0.0, 1.0);
                    double r = dis(gen);
                    if (r<var_freqs_rem[index2]) {
                        b.push_back(1);
                    } else {
                        b.push_back(0);
                    }
                    index2++;
                }
            }
            bitstrings.push_back(b);
        }
    } else {
        bitstrings=bitstrings_cut;
    }

}

vector< vector<int> > RandomString (int dim, int seq_length, int n, double& tolerance, const VectorXd& mu, const MatrixXd& omega, const VectorXi& positions) {
    SelfAdjointEigenSolver<MatrixXd> solver(omega);
    VectorXd eigenvalues = solver.eigenvalues();
    MatrixXd eigenvectors = solver.eigenvectors();
    for (int i=0;i<dim;i++) {
        if (eigenvalues(i) < tolerance) {
            eigenvalues(i)=0.0;
        }
    }
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> norm(0, 1);
    vector< vector<int> > random_strings;
    for (int s=0;s<n;s++) {
        VectorXd z(dim);
        for (int i=0;i<dim;i++) {
            z(i)=norm(gen);
        }
        VectorXd samplevector = mu + eigenvectors * eigenvalues.cwiseSqrt().asDiagonal() * z;
        vector<int> bitstring;
        for (int i=0;i<positions.size();i++) {
            double threshold = gsl_cdf_ugaussian_Pinv(1.0 - mu(i));
            if (samplevector(i)>threshold) {
                bitstring.push_back(1);
            } else {
                bitstring.push_back(0);
            }
        }
        random_strings.push_back(bitstring);
    }
    return random_strings;
}

double KingmanRatio (int n) {
    //Proportion more time in new branches
    if (n>1) {
        double ks=0;
        double ks2=0;
        for (int k=2;k<=n;k++) {
            ks=ks+(2./(k-1));
        }
        ks2=ks+(2./n);
        //cout << "Here " << ks << " " << ks2 << "\n";
        double prop=(ks2-ks)/ks;
        return prop;
    } else {
        return 0;
    }
}

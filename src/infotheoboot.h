#include <iostream>
//using namespace std;
#include <algorithm>
#include <ctime>
#include <cmath>
#include <vector>
#include <map>
#include <R.h>
#include <Rdefines.h>
#include <R_ext/Rdynload.h>

/* These functions are defined in the entropy.cpp file */

double digamma(double z);
double entropy_empirical(std::map< std::vector<int>,int > frequencies, int nb_samples);
double entropy_dirichlet(std::map<  std::vector<int> ,int > frequencies, int nb_samples, double beta);
double entropy_miller_madow(std::map<  std::vector<int> ,int > frequencies, int nb_samples);
double entropy_shrink(std::map<  std::vector<int> ,int > frequencies, int nb_samples);
double entropy(const int *d, int nsamples, int nvars, int c, bool *v);
double multiinformation(const int *d, int nsamples, int nvars, int c);
double interaction(const int *d, int nsamples, int nvars, int c);
double * megaentropy(const int *x, const int *y, const int *s,
  int perms, int nsamples, int xvars, int yvars, int c);
double entropy2(const int *dx, const int *dy, const int *ds,
  int nsamples, int i, int j, int c);
double entropyP(const int *dx, const int *dy, const int *ds,
    int nsamples, int i, int j, int c);

/* Entry points called from the R functions */
extern "C"
{
SEXP discEF( SEXP data, SEXP nrows, SEXP ncols, SEXP nbins );
SEXP discEW( SEXP data, SEXP nrows, SEXP ncols, SEXP nbins );
SEXP entropyR(SEXP data, SEXP nrows, SEXP ncols, SEXP choice);
SEXP bpcmiR(SEXP rxdata, SEXP rydata, SEXP rsdata, SEXP Rxcols, SEXP Rycols, SEXP Rrows, SEXP rpermreps, SEXP rchoi);
SEXP buildMIM(SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice);
SEXP multiinformationR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice);
SEXP interactionR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice);
}

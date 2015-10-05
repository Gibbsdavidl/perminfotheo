#include "infotheoboot.h"
#include <random>
#include <algorithm>
#include <iterator>

SEXP entropyR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
      SEXP res;
      PROTECT(Rdata = AS_INTEGER(Rdata));
      PROTECT(Rnrows= AS_INTEGER(Rnrows));
      PROTECT(Rncols= AS_INTEGER(Rncols));
	    PROTECT(Rchoice= AS_INTEGER(Rchoice));
      data = INTEGER_POINTER(Rdata);
      nrows= INTEGER_POINTER(Rnrows);
      ncols= INTEGER_POINTER(Rncols);
	    choice= INTEGER_POINTER(Rchoice);
      PROTECT(res = NEW_NUMERIC(1));
		  bool *sel = new bool[*ncols];
		  for( int i=0; i<*ncols; ++i )
			    sel[i] = true;
		  REAL(res)[0] = entropy(data, *nrows, *ncols, *choice, sel);
		  UNPROTECT(5);
      return res;
}

//res <- .Call( "bpcmiR",X,Y,S,n,m,r,permreps,choi,PACKAGE="infotheoboot")
SEXP bpcmiR(SEXP rxdata, SEXP rydata, SEXP rsdata,
  SEXP Rxcols, SEXP Rycols, SEXP Rrows,
  SEXP rpermreps, SEXP rchoi)
{
  const int *xdata, *ydata, *sdata;
  const int *nrows, *xcols, *ycols, *choice, *permreps;
  double *tmpres;
  SEXP res = PROTECT(allocVector(REALSXP, 3));
  PROTECT(rxdata= AS_INTEGER(rxdata));
  PROTECT(rydata= AS_INTEGER(rydata));
  PROTECT(rsdata= AS_INTEGER(rsdata)); // AS_Numeric
  PROTECT(Rrows= AS_INTEGER(Rrows));
  PROTECT(Rxcols= AS_INTEGER(Rxcols));
  PROTECT(Rycols= AS_INTEGER(Rycols));
  PROTECT(rpermreps= AS_INTEGER(rpermreps));
  PROTECT(rchoi= AS_INTEGER(rchoi));
  xdata = INTEGER_POINTER(rxdata);
  ydata = INTEGER_POINTER(rydata);
  sdata = INTEGER_POINTER(rsdata);
  nrows= INTEGER_POINTER(Rrows);
  xcols= INTEGER_POINTER(Rxcols);
  ycols= INTEGER_POINTER(Rycols);
  permreps= INTEGER_POINTER(rpermreps);
  choice= INTEGER_POINTER(rchoi);

  tmpres = megaentropy(xdata, ydata, sdata, *permreps, *nrows, *xcols, *ycols, *choice);

  REAL(res)[0] = tmpres[0];
  REAL(res)[1] = tmpres[1];
  REAL(res)[2] = tmpres[2];

  // free tmp res
  UNPROTECT(9);
  return res;
}


double * megaentropy(const int *x, const int *y, const int *s,
  int perms, int nsamples, int xvars, int yvars, int c) {

    double * res = new double [3];
    double Hsxy = 0.0, Hsx = 0.0, Hsy = 0.0, Hs = 0.0, obsSum = 0.0, permSum = 0.0, avgPerm = 0.0, pcount = 0.0;

    for (int i = 0; i < xvars; ++i) {
      for (int j = 0; j < yvars; ++j) {
        Hsxy = entropy2(x, y, s, nsamples, i,  j, c);
        Hsx  = entropy2(x, y, s, nsamples, i, -1, c);
        Hsy  = entropy2(x, y, s, nsamples, -1, j, c);
        Hs   = entropy2(x, y, s, nsamples, -1, -1, c);
        // sum entropy over all pairs
        obsSum += Hsy - Hs - Hsxy + Hsx;
      }
    }

    // for number of perms, and each pair therein
    for (int p = 0; p < perms; ++p) {
      permSum = 0.0;
      for (int i = 0; i < xvars; ++i) {
        for (int j = 0; j < yvars; ++j) {
          Hsxy = entropyP(x, y, s, nsamples, i,  j, c);
          Hsx  = entropyP(x, y, s, nsamples, i, -1, c);
          Hsy  = entropyP(x, y, s, nsamples, -1, j, c);
          Hs   = entropyP(x, y, s, nsamples, -1, -1, c);
          // sum entropy over all pairs
          permSum += Hsy - Hs - Hsxy + Hsx;
        }
      }
      // if sum greater than observed sum, increment perm count
      if (permSum > obsSum) {
        ++pcount;
      }
      avgPerm += permSum;
    }
    // return observed MI, avg perm MI, frac. perm greater
    res[0] = obsSum;
    res[1] = avgPerm/(double)perms;
    res[2] = pcount/(double)perms;
    return(res);
}


SEXP multiinformationR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
      SEXP res;
      PROTECT(Rdata = AS_INTEGER(Rdata));
      PROTECT(Rnrows= AS_INTEGER(Rnrows));
      PROTECT(Rncols= AS_INTEGER(Rncols));
		  PROTECT(Rchoice= AS_INTEGER(Rchoice));
      data = INTEGER_POINTER(Rdata);
      nrows= INTEGER_POINTER(Rnrows);
      ncols= INTEGER_POINTER(Rncols);
		  choice= INTEGER_POINTER(Rchoice);
      PROTECT(res = NEW_NUMERIC(1));
		  REAL(res)[0] = multiinformation(data, *nrows, *ncols, *choice);
		  UNPROTECT(5);
      return res;
}

SEXP interactionR (SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
         SEXP res;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows= AS_INTEGER(Rnrows));
         PROTECT(Rncols= AS_INTEGER(Rncols));
		 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);
		 choice= INTEGER_POINTER(Rchoice);
         PROTECT(res = NEW_NUMERIC(1));
		 REAL(res)[0] = interaction(data, *nrows, *ncols, *choice);
		 UNPROTECT(5);
      return res;
}

SEXP buildMIM(SEXP Rdata, SEXP Rnrows, SEXP Rncols, SEXP Rchoice)
{
      const int *data;
      const int *nrows, *ncols, *choice;
      double *ent, *res;
	  bool *sel;
         SEXP Rres;
         PROTECT(Rdata = AS_INTEGER(Rdata));
         PROTECT(Rnrows = AS_INTEGER(Rnrows));
         PROTECT(Rncols = AS_INTEGER(Rncols));
		 PROTECT(Rchoice= AS_INTEGER(Rchoice));
         data = INTEGER_POINTER(Rdata);
         nrows= INTEGER_POINTER(Rnrows);
         ncols= INTEGER_POINTER(Rncols);
		 choice= INTEGER_POINTER(Rchoice);
         PROTECT(Rres = NEW_NUMERIC((*ncols)*(*ncols)));
         res = NUMERIC_POINTER(Rres);
		 ent = new double[*ncols];
		 sel = new bool[*ncols];
		 for( int i=0; i<*ncols; ++i ){
			res[i*(*ncols)+i]=0;
			sel[i] = false;
		 }
		 for( int i=0; i<*ncols; ++i ){
			sel[i] = true;
			res[i*(*ncols)+i] = entropy(data, *nrows, *ncols, *choice, sel);
			sel[i] = false;
		 }
	     for( int i=1; i<*ncols; ++i ){
			sel[i] = true;
			for( int j=0; j<i; ++j ) {
				  sel[j] = true;
                  res[j*(*ncols)+i] = res[i*(*ncols)+j] = res[i*(*ncols)+i] + res[j*(*ncols)+j] - entropy(data, *nrows, *ncols, *choice, sel);
				  sel[j] = false;
			}
			sel[i] = false;
         }
         UNPROTECT(5);
      return Rres;
}


double digamma(double z) {
      if(z<=0) return 0;
      double zp, zpr, zprs, digam = 0;
         zp = z;
      while(zp < 30) {
             zpr = 1/zp;
             digam -= zpr;
             zp++;
      }
         zpr = 1/zp;
         zprs = zpr * zpr;
         digam += log(zp)+zpr*(-0.5+zpr*(-1.0/12.0+zprs*(1.0/120.0-zprs/252.0)));
      return digam;
}

double entropy_empirical(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e -= iter->second * log((double)iter->second);
      return log((double)nb_samples) + e/nb_samples;
}

double entropy_miller_madow(std::map< std::vector<int> ,int > frequencies, int nb_samples) {
      return entropy_empirical(frequencies,nb_samples) + (int(frequencies.size())-1)/(2.0*nb_samples);
}

double entropy_dirichlet(std::map< std::vector<int> ,int > frequencies, int nb_samples, double beta) {
      double e = 0;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            e+=(iter->second+beta)*(digamma(nb_samples+(frequencies.size()*beta)+1)-digamma(iter->second+beta+1));
      return e/(nb_samples+(frequencies.size()*beta));
}

double entropy_shrink(std::map< std::vector<int> ,int > frequencies, int nb_samples)
{
      double w = 0;
      int p = frequencies.size(), n2 = nb_samples*nb_samples;
      double lambda, beta;
      for (std::map< std::vector<int> ,int>::const_iterator iter = frequencies.begin(); iter != frequencies.end(); ++iter)
            w += iter->second*iter->second;
         lambda = p*(n2 - w)/((nb_samples-1)*(w*p - n2));
      if(lambda >= 1)
        return -log(1.0/p);
      else {
            beta = (lambda/(1-lambda))*nb_samples/frequencies.size();
        return entropy_dirichlet(frequencies, nb_samples, beta);
      }
}

double entropy2(const int *dx, const int *dy, const int *ds,
  int nsamples, int i, int j, int c) {

    std::map< std::vector<int> ,int > freq;
    std::vector<int> sel;
	  bool ok = true;
	  int nsamples_ok = 0;
    double H = 0;

    // for each sample
	  for(int s = 0; s < nsamples; ++s) {
		    ok = true;
		    sel.clear();

        if (i > -1) {  // take x?
          if(dx[s+i*nsamples]!=NA_INTEGER) {
               sel.push_back(dx[s+i*nsamples]);
          } else {
               ok = false;
          }
        }
        if (j > -1) {  // take y?
          if(dy[s+j*nsamples]!=NA_INTEGER) {
               sel.push_back(dy[s+j*nsamples]);
          } else {
               ok = false;
          }
        }
        if(ds[s]!=NA_INTEGER) {
             sel.push_back(ds[s]);
        } else {
             ok = false;
        }
		    if(ok) {
			       freq[sel]++;
			       nsamples_ok++;
	      }
	  } // end samples loop

	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}


double entropyP(const int *dx, const int *dy, const int *ds,
  int nsamples, int i, int j, int c) {

    std::map< std::vector<int> ,int > freq;
    std::vector<int> sel;
	  bool ok = true;
	  int nsamples_ok = 0;
    double H = 0;

    std::vector<int> idx1(nsamples);
    std::vector<int> idx2(nsamples);
    std::vector<int> idx3(nsamples);

    // fill the array
    for(int z = 0; z < nsamples; ++z){
      idx1.push_back(z);
      idx2.push_back(z);
      idx3.push_back(z);
    }

    std::random_device rd;
    std::mt19937 g(rd());

    if (i > -1) {  // take x?
      std::shuffle(idx1.begin(), idx1.end(), g);
    }
    if (j > -1) {  // take y?
      std::shuffle(idx2.begin(), idx2.end(), g);
    }
    std::shuffle(idx3.begin(), idx3.end(), g);

    // for each sample
	  for(int s = 0; s < nsamples; ++s) {
		    ok = true;
		    sel.clear();

        if (i > -1) {  // take x?
          if(dx[idx1.at(s)+i*nsamples]!=NA_INTEGER) {
               sel.push_back(dx[idx1.at(s)+i*nsamples]);
          } else {
               ok = false;
          }
        }
        if (j > -1) {  // take y?
          if(dy[idx2.at(s)+j*nsamples]!=NA_INTEGER) {
               sel.push_back(dy[idx2.at(s)+j*nsamples]);
          } else {
               ok = false;
          }
        }
        if(ds[s]!=NA_INTEGER) {
             sel.push_back(ds[idx3.at(s)]);
        } else {
             ok = false;
        }
		    if(ok) {
			       freq[sel]++;
			       nsamples_ok++;
	      }
	  } // end samples loop

	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}

double entropy(const int *d, int nsamples, int nvars, int c, bool *v) {
// H(d) using estimator c
	std::map< std::vector<int> ,int > freq;
	std::vector<int> sel;
	bool ok = true;
	int nsamples_ok = 0;
	double H = 0;
	for(int s = 0; s < nsamples; ++s) {
		ok = true;
		sel.clear();
		for(int i = 0; i < nvars; ++i) {
			if(v[i]) {
				if(d[s+i*nsamples]!=NA_INTEGER)
					sel.push_back(d[s+i*nsamples]);
				else
					ok = false;
			}
		}
		if(ok) {
			freq[sel]++;
			nsamples_ok++;
		}
	}
	if( c == 0 ) //empirical
		H = entropy_empirical(freq,nsamples_ok);
	else if( c == 1 ) //miller-madow
		H = entropy_miller_madow(freq,nsamples_ok);
	else if( c == 2 ) //dirichlet Schurmann-Grassberger
		H = entropy_dirichlet(freq,nsamples_ok, 1/freq.size());
	else if( c == 3 ) // shrink
		H = entropy_shrink(freq,nsamples_ok);
	return H;
}

double multiinformation(const int *d, int nsamples, int nvars, int c) {
		 bool *sel = new bool[nvars];
		 double sum = 0;
		 for( int i=0; i<nvars; ++i )
			sel[i] = false;
		 for(int i=0;i<nvars; ++i) {
			sel[i] = true;
			sum += entropy(d, nsamples, nvars, c, sel);
			sel[i] = false;
		 }
		 for( int i=0; i<nvars; ++i )
			sel[i] = true;
		 sum -= entropy(d, nsamples, nvars, c, sel);
		 return sum;
}



double interaction(const int *d, int nsamples, int nvars, int c) {
	double sum = 0; //
	int sign = 1;//
	bool *sel= new bool[nvars];//
	for(int i = 0; i < nvars; ++i)//
		sel[i] = false; //

	std::vector<int> indx;
	int n = nvars;
	int j=1;
	int k=n;
	for(int twk=j;twk<=k;twk++){
		int r=twk;
		bool done=true;
		for(int iwk=0;iwk<r;iwk++)indx.push_back(iwk);
		while(done){
			done=false;
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=true;	//
			sum += sign*entropy(d, nsamples, nvars, c, sel); //
			for(int owk=0;owk<r;owk++) //
				sel[indx[owk]]=false;	//
			for(int iwk=r-1;iwk>=0;iwk--){
				if(indx[iwk]<=(n-1)-(r-iwk)){
					indx[iwk]++;
					for(int swk=iwk+1;swk<r;swk++)
						indx[swk]=indx[swk-1]+1;
					iwk=-1;
					done=true;
				}
			}
		}
		sign = -sign; //
		indx.clear();
	}
	return sum;
}

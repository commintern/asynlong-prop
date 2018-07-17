#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;

/* ===========================================
 * Estimation for observation process 
 * ===========================================
 */

/*
 * The first part of the estimating equation of gamma
 */
// [[Rcpp::export]]
arma::vec ugamma1_C(Rcpp::ListOf<NumericMatrix>& kerMat,
                    Rcpp::ListOf<NumericMatrix>& covariates,
                    const unsigned int& n, const int& p) {
  unsigned int i = 0;
  arma::mat temp_kermat;
  arma::mat temp_cov,temp;
  arma::mat res = arma::zeros<arma::vec>(p);
  for (i = 0; i < n; i++) {
    temp_cov = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    temp = temp_kermat * temp_cov.t();
    res = res + arma::conv_to<arma::vec>::from(
      sum(temp, 0));
  }
  return res;
}

// Note that the definition of S^{(k)} is different from the paper here. For
// simiplicity, the divisor n is omitted
// Obtain Zbar=S^(1)/S^(1) for each T_{ij}. 
// [[Rcpp::export]]
Rcpp::List zbar_c(const arma::rowvec& gamma,
                  Rcpp::ListOf<NumericMatrix>& kerMat,
                  Rcpp::ListOf<NumericVector>& meas_times,
                  Rcpp::ListOf<NumericMatrix>& covariates,
                  const arma::vec& censor, const unsigned int& n,
                  const unsigned int& p) {
  unsigned int i, l;
  arma::dmat temp_kermat; 
  arma::mat temp, temp0, s0sumres, s1sumres;
  arma::mat tempres;
  arma::dmat temp_cov;
  arma::rowvec expgammaZ;
  arma::vec temp_meas_time;
  arma::vec censorind;

  List res(n);
  for (i = 0; i < n; i++) {
    temp_meas_time = vec(meas_times[i].begin(),meas_times[i].size(),false);
    s1sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, p);
    s0sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    for (l = 0; l < n; l++) {
      temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
      temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
      expgammaZ = exp(gamma * temp_cov);
      censorind = conv_to<arma::vec>::from(censor[l] > temp_meas_time);
      temp = (temp_kermat * (temp_cov.each_row() % expgammaZ).t());
      temp = temp.each_col() % censorind; 
      temp0 = temp_kermat * expgammaZ.t();
      temp0 = temp0.each_col() % censorind;
      //weightsum = sum(temp_kermat, 1);
      s1sumres = s1sumres + temp;
      s0sumres = s0sumres + temp0;
    }

    tempres = s1sumres.each_col() / s0sumres;
    tempres.transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
    res(i) = tempres.t();
  }
  return res;
}


// Second part of ugamma
// [[Rcpp::export]]
arma::vec ugamma2_C(const arma::rowvec& gamma,
                    Rcpp::ListOf<NumericMatrix>& kerMat,
                    Rcpp::ListOf<NumericVector>& meas_times,
                    Rcpp::ListOf<NumericMatrix>& covariates,
                    const arma::vec& censor, const unsigned int& n,
                    const unsigned int& p) {
  unsigned int i = 0;
  List zbarres = zbar_c(gamma, kerMat, meas_times, covariates, censor, n, p);
  arma::mat temp_kermat;
  arma::mat temp_zbari;
  arma::mat res = arma::zeros<arma::vec>(p);
  for (i = 0; i < n; i++) {
    temp_zbari = Rcpp::as<arma::mat>(zbarres[i]);
    //mat(zbarres[i].begin(),zbarres[i].nrow(),zbarres[i].ncol(),false);
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    
    res = res + arma::conv_to<arma::vec>::from(sum(temp_zbari * temp_kermat, 1));
  }
  return res;
}


// Calculate S^{(0)} for \Lambda(t)
// TODO change the function name this is not S0_C
// [[Rcpp::export]]
Rcpp::List dlambda_C(const arma::rowvec& gamma,
                Rcpp::ListOf<NumericMatrix>& kerMat,
                Rcpp::ListOf<NumericVector>& meas_times,
                Rcpp::ListOf<NumericMatrix>& covariates,
                const arma::vec& censor, const unsigned int& n,
                const unsigned int& p) {
  unsigned int i, l;
  arma::dmat temp_kermat, temp, temp0, s0sumres;
  arma::dmat temp_cov;
  arma::mat tempres;
  arma::rowvec expgammaZ;
  arma::vec temp_meas_time;
  arma::vec censorind;
  List res(n);
  for (i = 0; i < n; i++) {
    temp_meas_time = Rcpp::as<arma::vec>(meas_times[i]);
    s0sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    //s0weightsum = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    for (l = 0; l < n; l++) {
      temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
      temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
      expgammaZ = exp(gamma * temp_cov);
      censorind = conv_to<arma::vec>::from(censor[l] > temp_meas_time);
      temp0 =  temp_kermat * expgammaZ.t();
      temp0 = temp0.each_col() % censorind;
      //s0weightsum = s0weightsum + sum(temp_kermat, 1);
      s0sumres = s0sumres + temp0;
    }
    // s1(i) = s1sumres;
    // s0(i) = s0sumres;
    tempres = sum(mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false),1)/s0sumres;
    tempres.transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
    res(i) =  tempres;
  }
  return res;
}




/*** R
### Test
#ugamma_C(simdatasam,0.5,log)
#foo(kerMat)
#ugamma1_C(kerMat,covarlist,2L,1)
#zbar(1,kerMat,timelist,covarlist,censorlist,2,2)
#zbar(1,kerMat,timelist,covarlist,censorlist,2,2)

*/

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//#include "optim.hpp"
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;
#include "numeric.hpp"
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
  arma::mat res = arma::zeros<arma::mat>(p,1);
  for (i = 0; i < n; i++) {
    temp_cov = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    res = res + sum(temp_kermat * temp_cov.t(), 0);
  }
  return conv_to<arma::vec>::from(res);
}

// Note that the definition of S^{(k)} is different from the paper here. For
// simiplicity, the divisor n is omitted
// Obtain Zbar=S^(1)/S^(1) for each T_{ij}.
// [[Rcpp::export]]
arma::field<arma::mat> zbar_c(const arma::rowvec& gamma,
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

  arma::field<arma::mat>  res(n);
  for (i = 0; i < n; i++) {
    temp_meas_time = vec(meas_times[i].begin(),meas_times[i].size(),false);
    s1sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, p);
    s0sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    for (l = 0; l < n; l++) {
      temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
      temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
      expgammaZ = exp(gamma * temp_cov);
      censorind = conv_to<arma::vec>::from(censor[l] > temp_meas_time);
      cout << censorind << endl;
      temp = (temp_kermat * (temp_cov.each_row() % expgammaZ).t());
      temp = temp.each_col() % censorind;
      temp0 = temp_kermat * expgammaZ.t();
      temp0 = temp0.each_col() % censorind;
      //weightsum = sum(temp_kermat, 1);
      s1sumres = s1sumres + temp;
      s0sumres = s0sumres + temp0;
    }

    tempres = s1sumres.each_col() / s0sumres;
    tempres.replace(datum::nan, 0);
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
  unsigned int i,l = 0;
  //arma::field<arma::mat> zbarres = zbar_c(gamma, kerMat, meas_times, covariates, censor, n, p);
  arma::mat temp_kermat,temp_cov,temp_KerexpgamZ,temp_expgamZ;

  //arma::mat temp_zbari;
  //arma::field < arma::mat > zbar_list(n);
  //arma::mat temp;
  arma::vec temp_meas_time_i;
  arma::mat zbar_temp_i,temp;
  arma::vec censorind;
  arma::rowvec expgammaZ;
  arma::vec res = arma::zeros<arma::vec>(p);
  arma::vec den_temp;
  for (i = 0; i < n; i++) {
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    zbar_temp_i = arma::zeros<arma::dmat>(temp_meas_time_i.n_elem, p);
    den_temp = arma::zeros<arma::vec>(temp_meas_time_i.n_elem);
    for( l =0; l<n; l++){
      //temp_zbari = Rcpp::as<arma::mat>(zbarres[i]);
      //mat(zbarres[i].begin(),zbarres[i].nrow(),zbarres[i].ncol(),false);
      temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
      temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
      temp_expgamZ = exp(gamma * temp_cov);
      //temp_KerexpgamZ = temp_kermat.each_row() % exp(gamma * temp_cov);
      censorind = conv_to < arma::vec > ::from(censor[l] > temp_meas_time_i);
      //temp_KerexpgamZ.each_col() %= censorind;
      temp = temp_kermat * (temp_cov.each_row() % temp_expgamZ).t();
      zbar_temp_i += temp.each_col() %  censorind;

      den_temp += temp_kermat * temp_expgamZ.t() % censorind;
    }
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    zbar_temp_i.each_col() /= den_temp;
    res += zbar_temp_i.t() * sum(temp_kermat,1);
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
    //temp_meas_time = Rcpp::as<arma::vec>(meas_times[i]);
    temp_meas_time = vec(meas_times[i].begin(),meas_times[i].size(),false);
    s0sumres = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    //s0weightsum = arma::zeros<arma::dmat>(temp_meas_time.n_elem, 1);
    for (l = 0; l < n; l++) {
      temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
      temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
      //expgammaZ = exp(gamma * temp_cov);
      censorind = conv_to<arma::vec>::from(censor[l] > temp_meas_time);
      temp0 =  temp_kermat * exp(gamma * temp_cov).t();
      temp0.each_col() %= censorind;
      //s0weightsum = s0weightsum + sum(temp_kermat, 1);
      s0sumres = s0sumres + temp0;
    }
    // s1(i) = s1sumres;
    // s0(i) = s0sumres;
    tempres = sum(mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false),1)/s0sumres;
    tempres.replace(datum::nan, 0);
    res(i) =  tempres;
  }
  return res;
}

/* ==================================================================================================
 * ===================================================================================================
 */

// [[Rcpp::export]]
arma::vec ugamma_C_1(const arma::rowvec& gamma,
                     Rcpp::ListOf<NumericMatrix>& kerMat,
                     Rcpp::ListOf<NumericVector>& meas_times,
                     Rcpp::ListOf<NumericMatrix>& covariates,
                     const arma::vec& censor, const unsigned int& n,
                     const unsigned int& p) {
  unsigned int i,j,k,l,kk,Jn,Kn, KKn = 0;
  arma::vec temp_meas_time;
  arma::mat temp_kermat,temp_kermat_il;
  arma::mat temp_cov,temp,temp_cov_l;
  arma::mat res = arma::zeros<arma::mat>(p,1);

  double temp_w;

  arma::vec temp_S1;
  double temp_S0;
  arma::vec temp_Zbar;
  for (i = 0; i < n; i++) {
    temp_cov = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    temp_meas_time = vec(meas_times[i].begin(), meas_times[i].size(), false);
    Jn = temp_kermat.n_cols;
    Kn = temp_kermat.n_rows;


    // Obtain gamma
    for(j = 0;j<Jn; j++){

      // Obtain Z_bar
      temp_S1 = arma::zeros<arma::vec>(p);
      temp_S0 = 0;
      for(l = 0; l <n; l++){
        if(censor[l] > temp_meas_time[j]){
          temp_kermat_il = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
          temp_cov_l = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
          KKn = temp_kermat_il.n_cols;
          for(kk = 0; kk <KKn; kk++){
            if(abs(temp_kermat_il.at(j,kk)-1e-10)>0){
              temp_w=  temp_kermat_il.at(j,kk)*exp(as_scalar(gamma * temp_cov_l.col(kk)));
              temp_S1 = temp_S1  + temp_w*temp_cov_l.col(kk);
              temp_S0 = temp_S0  + temp_w;
            }
          }
        }
      }

      temp_Zbar = temp_S1/temp_S0;

      // Summation
      for(k = 0;j<Kn; k++){
        if(abs(temp_kermat.at(j,k)-1e-10)>0){
          res = res + temp_kermat.at(j,k)*(temp_cov.col(k)-temp_Zbar);
        }
      }
    }


  }
  return conv_to<arma::vec>::from(res);
}



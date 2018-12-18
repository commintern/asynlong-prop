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
// // [[Rcpp::export]]
// arma::vec ugamma2_C(const arma::rowvec& gamma,
//                     Rcpp::ListOf<NumericMatrix>& kerMat,
//                     Rcpp::ListOf<NumericVector>& meas_times,
//                     Rcpp::ListOf<NumericMatrix>& covariates,
//                     const arma::vec& censor, const unsigned int& n,
//                     const unsigned int& p) {
//   unsigned int i = 0;
//   arma::field<arma::mat> zbarres = zbar_c(gamma, kerMat, meas_times, covariates, censor, n, p);
//   arma::mat temp_kermat;
//   arma::mat temp_zbari;
//   arma::mat res = arma::zeros<arma::vec>(p);
//   for (i = 0; i < n; i++) {
//     //temp_zbari = Rcpp::as<arma::mat>(zbarres[i]);
//     //mat(zbarres[i].begin(),zbarres[i].nrow(),zbarres[i].ncol(),false);
//     temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
//
//     res = res + arma::conv_to<arma::vec>::from(zbarres[i] * sum(temp_kermat, 1));
//   }
//   return res;
// }


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
      zbar_temp_i += temp.each_col() /  censorind;

      den_temp += temp_kermat * temp_expgamZ.t() / censorind;
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


//==========================================================================================================

// struct sim_data_t
// {
//   Rcpp::ListOf<NumericMatrix>* kerMat;
//   Rcpp::ListOf<NumericVector>* meas_times;
//   Rcpp::ListOf<NumericMatrix>* covariates;
//   arma::vec ugamma1;
//   const arma::vec* censor;
//   const unsigned int* n;
//   const unsigned int* p;
// };
//
// arma::vec ugamma_wrapper_C(const arma::vec& gamma_inp, void* opt_data){
//   sim_data_t* objfn_data = reinterpret_cast<sim_data_t*>(opt_data);
//   arma::rowvec gamma_inp_row = gamma_inp.t();
//   arma::vec ugamma2 = ugamma2_C(gamma_inp_row,
//                             *(objfn_data->kerMat),
//                             *(objfn_data->meas_times),
//                             *(objfn_data->covariates),
//                             *(objfn_data->censor),
//                             *(objfn_data->n),
//                             *(objfn_data->p));
//   return objfn_data->ugamma1-ugamma2;
// }
//
// struct Ugamma_fun
// {
//   Rcpp::ListOf<NumericMatrix>* kerMat;
//   Rcpp::ListOf<NumericVector>* meas_times;
//   Rcpp::ListOf<NumericMatrix>* covariates;
//   const arma::vec* censor;
//   const unsigned int* n;
//   const unsigned int* p;
//   arma::vec ugamma1;
//   void cal_ugamma1() { ugamma1 =  ugamma1_C(*kerMat,*covariates,*n,*p);}
//   arma::vec operator()(arma::vec gamma) {
//     return ugamma1 - ugamma2_C(gamma, *kerMat, *meas_times, *covariates, *censor, *n, *p);
//     //return pow(x,3) -gamma;
//   }
// };

// // [[Rcpp::export]]
// arma::vec ugamma_wrapper_test_C(const arma::rowvec& gamma,
//                                  Rcpp::ListOf<NumericMatrix>& kerMat,
//                                  Rcpp::ListOf<NumericVector>& meas_times,
//                                  Rcpp::ListOf<NumericMatrix>& covariates,
//                                  const arma::vec& censor, const unsigned int& n,
//                                  const unsigned int& p){
//   sim_data_t objfn_data;
//   objfn_data.kerMat = &kerMat;
//   objfn_data.meas_times = &meas_times;
//   objfn_data.covariates = &covariates;
//   objfn_data.censor = &censor;
//   objfn_data.n = &n;
//   objfn_data.p = &p;
//   //arma::rowvec gamma_inp = gamma;
//   arma::vec res = ugamma_wrapper_C(gamma, &objfn_data);
//   return res;
// }

//
// // [[Rcpp::export]]
// arma::vec gammaest_C(arma::vec& gamma_inp,
//                      Rcpp::ListOf<NumericMatrix>& kerMat,
//                      Rcpp::ListOf<NumericVector>& meas_times,
//                      Rcpp::ListOf<NumericMatrix>& covariates,
//                      const arma::vec& censor, const unsigned int& n,
//                      const unsigned int& p){
//   arma::vec ugamma1 = ugamma1_C(kerMat,covariates,n,p);
//   sim_data_t objfn_data;
//   objfn_data.kerMat = &kerMat;
//   objfn_data.ugamma1 = ugamma1;
//   objfn_data.meas_times = &meas_times;
//   objfn_data.covariates = &covariates;
//   objfn_data.censor = &censor;
//   objfn_data.n = &n;
//   objfn_data.p = &p;
//   arma::vec gamma_inp_temp = gamma_inp;
//   bool success = optim::broyden_df(gamma_inp_temp, ugamma_wrapper_C, &objfn_data);
//   return gamma_inp_temp;
// }

//// [[Rcpp::export]]
// arma::vec gammaest_C(arma::vec& gamma_inp,
//                      Rcpp::ListOf<NumericMatrix>& kerMat,
//                      Rcpp::ListOf<NumericVector>& meas_times,
//                      Rcpp::ListOf<NumericMatrix>& covariates,
//                      const arma::vec& censor, const unsigned int& n,
//                      const unsigned int& p, arma::vec& ugamma1){
//   //arma::vec ugamma1 = ugamma1_C(kerMat,covariates,n,p);
//   Ugamma_fun ugamma;
//   ugamma.kerMat = &kerMat;
//   ugamma.ugamma1 = ugamma1;
//   ugamma.meas_times = &meas_times;
//   ugamma.covariates = &covariates;
//   ugamma.censor = &censor;
//   ugamma.n = &n;
//   ugamma.p = &p;
//   //ugamma.cal_ugamma1();
//   arma::vec gamma_inp_temp = gamma_inp;
//   bool check = false;
//   arma::vec res = broyden(gamma_inp_temp, ugamma, check);
//   return res;
// }

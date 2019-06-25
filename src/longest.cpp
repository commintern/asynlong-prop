#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <iostream>


using namespace arma;
using namespace Rcpp;
using namespace std;



/* ===========================================
 * Estimation for longitudinal processes with additive modelgpopo
 * ===========================================
 */

// [[Rcpp::export]]
arma::vec countprofun_C(const arma::vec& counttime, const arma::vec& externalTime){
  arma::uword i,j;
  arma::uword Nc = counttime.size();
  arma::uword Nt = externalTime.size();
  arma::vec res(Nt);
  j = 0;
  i = 0;
  while((i<Nt)&&(j<Nc)){
    if(j>=Nc-1){
      if(counttime(j)<externalTime(i)){
        res(i) = Nc;
      } else {
        res(i) = Nc-1;
      }
      i++;
    } else {
      if(counttime(j)<externalTime(i)){
        if(counttime(j+1)>externalTime(i)){
          res(i) = j+1;
          ++i;
        } else {
          ++j;
        }
      } else {
        res(i) = j;
        ++i;
      }
    }

  }
  return res;
}

// [[Rcpp::export]]
arma::cube Xgen_C(const arma::mat& covMat, const arma::vec& countprocess, const unsigned int& p){
  arma::uword i;
  arma::uword nrow = countprocess.size();
  arma::uword ncol = covMat.n_cols;
  arma::cube res(nrow,ncol,p+1);
  // for(i=0;i<nrow;i++){
  //   res.subcube(i,0,0,i,ncol-1,p-1) = covMat;
  // }
  for(i=0;i<p;i++){
    res.slice(i).each_row() = covMat.row(i);
  }
  res.slice(p).each_col() = countprocess;

  return res;
}



// [[Rcpp::export]]
arma::vec longest_prop_c(const arma::rowvec & theta,
                         const arma::rowvec & gamma,
                     Rcpp::ListOf < NumericMatrix > & kerMat,
                     Rcpp::ListOf < NumericVector > & meas_times,
                     Rcpp::ListOf < NumericMatrix > & covariates,
                     Rcpp::ListOf < NumericVector > & response,
                     Rcpp::ListOf < NumericVector > & dlambda,
                     const arma::vec & censor,
                     const unsigned int & n,
                     const unsigned int & p) {
  unsigned int i, l = 0;
  arma::mat temp_kermat;

  arma::rowvec expgammaZ;
  arma::field < arma::rowvec > expgamZ_list(n);
  arma::field < arma::mat > Xbar_list(n);
  arma::field < arma::mat > KerexpgamZ_list(n * n);
  arma::field < arma::cube > Xmat_list(n * n), KerXexpgamZ_list(n * n),  XXtbar_list(n);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;
  arma::vec censorind, temp_response;

  unsigned int j, k, Jn, Kn;
  // Rcpp::Rcout <<  "start" <<endl;
  // Calculate Xbar and XXtbar
  // arma::cube Xmat_temp;
  arma::mat temp_cov, temp_KerexpgamZ;
  arma::rowvec tempexpgamZ;
  // arma::mat  XXtbar_num_temp=mat(p+1,p+1,arma::fill::zeros);
  arma::mat temp_kermat_ii;
  // arma::vec  Xbar_num_temp=vec(p+1,arma::fill::zeros);
  arma::vec den_temp;
  //Rcpp::Rcout <<  "00" <<endl;
  arma::vec temp_X = vec(p + 1, arma::fill::zeros);
  arma::vec thetanum_part1 = vec(p + 1, arma::fill::zeros);
  arma::vec temp_thetanum_part1 = vec(p + 1, arma::fill::zeros);
  arma::vec thetanum_part2 = vec(p + 1, arma::fill::zeros);
  //arma::mat thetaden = mat(p + 1, p + 1, arma::fill::zeros);
  arma::vec temp_kermat_rowsum;
  // Rcpp::Rcout <<  "11" <<endl;
  // Calculate \hat{theta}
  for (i = 0; i < n; i++) {
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    temp_response = vec(response[i].begin(), response[i].size(), false);
    Jn = temp_meas_time_i.size();
    Xbar_list(i) = mat(Jn, p + 1, arma::fill::zeros);
    //XXtbar_list(i) = cube(p + 1, p + 1, Jn, arma::fill::zeros);
    // temp_thetanum_part1 = vec(p+1,arma::fill::zeros);
    den_temp = vec(Jn, arma::fill::zeros);
    // Rcpp::Rcout <<  "22" <<endl;
    for (l = 0; l < n; l++) {
      temp_meas_time_l =
        vec(meas_times[l].begin(), meas_times[l].size(), false);
      temp_cov = mat(covariates[l].begin(), covariates[l].nrow(),
                     covariates[l].ncol(), false);
      temp_kermat = mat(kerMat[i * n + l].begin(), kerMat[i * n + l].nrow(),
                        kerMat[i * n + l].ncol(), false);
      temp_KerexpgamZ = temp_kermat.each_row() % exp(gamma * temp_cov);
      ///////////////////////////////////// TODO
      censorind = conv_to < arma::vec > ::from(censor[l] > temp_meas_time_i);
      temp_KerexpgamZ.each_col() %= censorind;
      temp_countprocess_l = countprofun_C(temp_meas_time_l, temp_meas_time_i);

      //den_temp += sum(temp_KerexpgamZ, 1);

      Kn = temp_cov.n_cols;
      // Rcpp::Rcout <<  "33" <<endl;

      for (j = 0; j < Jn; j++) {
        temp_thetanum_part1.zeros();
        for (k = 0; k < Kn; k++) {
          temp_X.head(p) = temp_cov.col(k);
          temp_X(p) = temp_countprocess_l(j);
          Xbar_list(i).row(j) += temp_KerexpgamZ(j, k) * exp(theta * temp_X) * temp_X.t();
          den_temp.at(j) += as_scalar(temp_KerexpgamZ(j, k) * exp(theta * temp_X));
          //XXtbar_list(i).slice(j) +=
          //  temp_KerexpgamZ(j, k) * temp_X * temp_X.t();
          // Calculate the numerator of theta equation part1;
          if (i == l) {
            temp_thetanum_part1 += temp_kermat(j, k) * temp_X;
            // thetanum_part1+= temp_kermat(j,k)*temp_X*temp_response(j);
          }
        }
        if (i == l) {
          thetanum_part1 += temp_thetanum_part1 * temp_response(j);
        }
      }
    }

    // Rcpp::Rcout <<  "44" <<endl;
    Xbar_list(i).each_col() /= den_temp;
    //for (j = 0; j < Jn; j++) {
    //  XXtbar_list(i).slice(j) /= den_temp(j);
    //}

    // Rcpp::Rcout <<  "55" <<endl;
    temp_kermat_ii = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                         kerMat[i * n + i].ncol(), false);
    temp_kermat_rowsum = sum(temp_kermat_ii, 1);
    thetanum_part2 +=
      sum(Xbar_list(i).each_col() % (temp_response % temp_kermat_rowsum), 0).t();
      // for (j = 0; j < Jn; j++) {
      //   thetaden += (XXtbar_list(i).slice(j) -
      //     Xbar_list(i).row(j).t() * Xbar_list(i).row(j)) *
      //     temp_kermat_rowsum(j);
      // }
      // thetaden += sum(Xbar_list(i).each_col() % response(i) %
      // sum(temp_kermat_ii,1),0).t();
  }

  arma::vec thetanum = thetanum_part1 - thetanum_part2;
  //arma::vec thetaest = inv_sympd(thetaden) * thetanum;
  //arma::field < arma::vec > gmu0est(n);
  //for (i = 0; i < n; i++) {
  //  temp_response = vec(response[i].begin(), response[i].size(), false);
  //  gmu0est(i) = temp_response - Xbar_list(i) * thetaest;
  //}
  return thetanum;

  //   return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
  //                             Rcpp::Named("Xbar_list") = Xbar_list,
  //                             Rcpp::Named("temp0") = thetaden,
  //                             Rcpp::Named("temp1") = thetanum  );
}

// ============================ Purturbation ==================================================

// [[Rcpp::export]]
Rcpp::List longest_pur_c(const arma::rowvec & gamma,
                     Rcpp::ListOf < NumericMatrix > & kerMat,
                     Rcpp::ListOf < NumericVector > & meas_times,
                     Rcpp::ListOf < NumericMatrix > & covariates,
                     Rcpp::ListOf < NumericVector > & response,
                     Rcpp::ListOf < NumericVector > & dlambda,
                     const arma::vec & censor,
                     const unsigned int & n,
                     const unsigned int & p, const arma::vec& pur_weights) {
  unsigned int i, l = 0;
  arma::mat temp_kermat;

  arma::rowvec expgammaZ;
  arma::field < arma::rowvec > expgamZ_list(n);
  arma::field < arma::mat > Xbar_list(n);
  arma::field < arma::mat > KerexpgamZ_list(n * n);
  arma::field < arma::cube > Xmat_list(n * n), KerXexpgamZ_list(n * n),  XXtbar_list(n);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;
  arma::vec censorind, temp_response;

  unsigned int j, k, Jn, Kn;
  // Rcpp::Rcout <<  "start" <<endl;
  // Calculate Xbar and XXtbar
  // arma::cube Xmat_temp;
  arma::mat temp_cov, temp_KerexpgamZ;
  arma::rowvec tempexpgamZ;
  // arma::mat  XXtbar_num_temp=mat(p+1,p+1,arma::fill::zeros);
  arma::mat temp_kermat_ii;
  // arma::vec  Xbar_num_temp=vec(p+1,arma::fill::zeros);
  arma::vec den_temp;
  //Rcpp::Rcout <<  "00" <<endl;
  arma::vec temp_X = vec(p + 1, arma::fill::zeros);
  arma::vec thetanum_part1 = vec(p + 1, arma::fill::zeros);
  arma::vec temp_thetanum_part1 = vec(p + 1, arma::fill::zeros);
  arma::vec thetanum_part2 = vec(p + 1, arma::fill::zeros);
  arma::mat thetaden = mat(p + 1, p + 1, arma::fill::zeros);
  arma::vec temp_kermat_rowsum;
  // Rcpp::Rcout <<  "11" <<endl;
  // Calculate \hat{theta}
  for (i = 0; i < n; i++) {
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    temp_response = vec(response[i].begin(), response[i].size(), false);
    Jn = temp_meas_time_i.size();
    Xbar_list(i) = mat(Jn, p + 1, arma::fill::zeros);
    XXtbar_list(i) = cube(p + 1, p + 1, Jn, arma::fill::zeros);
    // temp_thetanum_part1 = vec(p+1,arma::fill::zeros);
    den_temp = vec(Jn, arma::fill::zeros);
    // Rcpp::Rcout <<  "22" <<endl;
    for (l = 0; l < n; l++) {
      temp_meas_time_l =
        vec(meas_times[l].begin(), meas_times[l].size(), false);
      temp_cov = mat(covariates[l].begin(), covariates[l].nrow(),
                     covariates[l].ncol(), false);
      temp_kermat = mat(kerMat[i * n + l].begin(), kerMat[i * n + l].nrow(),
                        kerMat[i * n + l].ncol(), false);
      temp_KerexpgamZ = temp_kermat.each_row() % exp(gamma * temp_cov);
      ///////////////////////////////////// TODO
      censorind = conv_to < arma::vec > ::from(censor[l] > temp_meas_time_i);
      temp_KerexpgamZ.each_col() %= censorind;
      temp_countprocess_l = countprofun_C(temp_meas_time_l, temp_meas_time_i);

      den_temp += sum(temp_KerexpgamZ, 1);

      Kn = temp_cov.n_cols;
      // Rcpp::Rcout <<  "33" <<endl;

      for (j = 0; j < Jn; j++) {
        temp_thetanum_part1.zeros();
        for (k = 0; k < Kn; k++) {
          temp_X.head(p) = temp_cov.col(k);
          temp_X(p) = temp_countprocess_l(j);
          Xbar_list(i).row(j) += temp_KerexpgamZ(j, k) * temp_X.t();
          XXtbar_list(i).slice(j) +=
            temp_KerexpgamZ(j, k) * temp_X * temp_X.t();
          // Calculate the numerator of theta equation part1;
          if (i == l) {
            temp_thetanum_part1 += temp_kermat(j, k) * temp_X;
            // thetanum_part1+= temp_kermat(j,k)*temp_X*temp_response(j);
          }
        }
        if (i == l) {
          thetanum_part1 += pur_weights(i)*temp_thetanum_part1 * temp_response(j);
        }
      }
    }

    // Rcpp::Rcout <<  "44" <<endl;
    Xbar_list(i).each_col() /= den_temp;
    for (j = 0; j < Jn; j++) {
      XXtbar_list(i).slice(j) /= den_temp(j);
    }

    // Rcpp::Rcout <<  "55" <<endl;
    temp_kermat_ii = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                         kerMat[i * n + i].ncol(), false);
    temp_kermat_rowsum = sum(temp_kermat_ii, 1);
    thetanum_part2 +=
      pur_weights(i)*sum(Xbar_list(i).each_col() % (temp_response % temp_kermat_rowsum), 0)
                      .t();
      for (j = 0; j < Jn; j++) {
        thetaden += pur_weights(i)*(XXtbar_list(i).slice(j) -
          Xbar_list(i).row(j).t() * Xbar_list(i).row(j)) *
          temp_kermat_rowsum(j);
      }
      // thetaden += sum(Xbar_list(i).each_col() % response(i) %
      // sum(temp_kermat_ii,1),0).t();
  }

  arma::vec thetanum = thetanum_part1 - thetanum_part2;
  arma::vec thetaest = inv_sympd(thetaden) * thetanum;
  arma::field < arma::vec > gmu0est(n);
  for (i = 0; i < n; i++) {
    temp_response = vec(response[i].begin(), response[i].size(), false);
    gmu0est(i) = temp_response - Xbar_list(i) * thetaest;
  }
  return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest);

  //   return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
  //                             Rcpp::Named("Xbar_list") = Xbar_list,
  //                             Rcpp::Named("temp0") = thetaden,
  //                             Rcpp::Named("temp1") = thetanum  );
}




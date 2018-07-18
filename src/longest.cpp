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
arma::cube Xgen_C(const arma::mat& covMat, const arma::vec& countprocess, const int& p){
  arma::uword i;
  arma::uword nrow = countprocess.size();
  arma::uword ncol = covMat.n_cols;
  arma::cube res(nrow,ncol,p+1);
  for(i=0;i<nrow;i++){
    res.subcube(i,0,0,i,ncol-1,p-1) = covMat;
  }
  res.slice(p).each_col() = countprocess;

  return res;
}

/*
 * \exp\left(\gamma^{\prime}Z_{i}\left(R_{ik}\right)\right)

 * K\left(T_{lu}-R_{ik}\right)I\left(C_{i}\ge T_{lu}\right)\exp\left(\gamma^{\prime}Z_{i}\left(R_{ik}\right)\right)

 * K\left(T_{lu}-R_{ik}\right)I\left(C_{i}\ge t\right)X_{i}\left(R_{ik},T_{lu}\right)\exp\left(\gamma^{\prime}Z_{i}\left(R_{ik}\right)\right)
 */

// [[Rcpp::export]]
Rcpp::List longest_c(const arma::rowvec& gamma,
                  Rcpp::ListOf<NumericMatrix>& kerMat,
                  Rcpp::ListOf<NumericVector>& meas_times,
                  Rcpp::ListOf<NumericMatrix>& covariates,
                  Rcpp::ListOf<NumericVector>& response,
                  Rcpp::ListOf<NumericVector>& dlambda,
                  const arma::vec& censor, const unsigned int& n,
                  const unsigned int& p) {
  unsigned int i, l=0;
  arma::mat temp_kermat,temp,temp_test;
  arma::mat temp1, temp0, s0sumres, s1sumres;
  arma::cube temp_cube;
  arma::mat tempres;
  double temp_dvalue;
  arma::dmat temp_cov_i;
  arma::rowvec expgammaZ;
  arma::field<arma::rowvec> expgamZ_list(n);
  arma::field<arma::mat> Xbar_list(n);
  arma::field<arma::mat> KerexpgamZ_list(n*n);
  arma::field<arma::cube> Xmat_list(n*n),KerXexpgamZ_list(n*n), XXtbar_list(n);
  arma::vec temp_meas_time_i,temp_meas_time_l,temp_countprocess_l,temp_vec1,temp_vec2;
  arma::vec censorind;

  unsigned int j,k,Jn,Kn;
  //Rcpp::Rcout <<  "start" <<endl;
  //List res(2);
  // i is the iterator for R_i
  for (i = 0; i < n; i++) {
    temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
    temp_cov_i = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
    // Obtain exp(gamma * Z_i(R_{ik}))
    // Store them in expgamZ_list

    expgamZ_list(i) = conv_to<arma::rowvec>::from(exp(gamma * temp_cov_i));

    // l is the iterator for T_l
    for (l = 0; l < n; l++) {
      // K(T_{lu}-R_{ik})I(C_i>=T_lu)*exp(gamma * Z_i(R_{ik}))
      // Store them in KerexpgamZ_list

      temp_kermat = mat(kerMat[l * n + i].begin(),kerMat[l * n + i].nrow(),kerMat[l * n + i].ncol(),false);
      temp_meas_time_l = vec(meas_times[l].begin(),meas_times[l].size(),false);
      KerexpgamZ_list(l*n + i) = temp_kermat.each_row() % expgamZ_list(i);
      censorind = conv_to<arma::vec>::from(censor[i] > temp_meas_time_l);
      KerexpgamZ_list(l*n + i).each_col() %= censorind;


      // Genrate N_i(T_{lu})
      temp_countprocess_l = countprofun_C(temp_meas_time_i, temp_meas_time_l);

      // Get X_i(R_(ik),T_(lu))
      Xmat_list(l*n + i)  = Xgen_C(temp_cov_i,temp_countprocess_l,p);

      // Calculate
      //K\left(T_{lu}-R_{ik}\right)I\left(C_{i}\ge t\right)X_{i}\left(R_{ik},T_{lu}\right)
      //  *\exp\left(\gamma^{\prime}Z_{i}\left(R_{ik}\right)\right)
      KerXexpgamZ_list(l*n+i) = Xmat_list(l*n + i).each_slice() % KerexpgamZ_list(l*n + i);
    }
  }

  //Rcpp::Rcout << expgamZ_list <<endl;
  //Rcpp::Rcout << KerexpgamZ_list <<endl;
  //Rcpp::Rcout <<  "2step" <<endl;

  // Calculate \bar{X}
  for(i=0;i<n;i++){
    temp1 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, p+1);
    temp0 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, 1);

    for(l=0;l<n;l++){
      temp = arma::sum(KerXexpgamZ_list(i*n + l),1);
      temp1 = temp1 + temp;
      temp0 = temp0 + arma::sum(KerexpgamZ_list(i*n + l),1);
    }
    Xbar_list(i) = temp1.each_col() /temp0;
    //Xbar_list(i).transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
  }

  // TODO futher explain
  // Calculate \bar{XXt}
  //Rcpp::Rcout <<  "2.5step" <<endl;
  for(i=0;i<n;i++){
    //temp_cube = arma::cube(p+1,p+1,KerexpgamZ_list(i*n+0).n_rows);
    Jn = KerexpgamZ_list(i*n +0).n_rows;
    XXtbar_list(i) = arma::cube(p+1,p+1,Jn);
    for(j=0;j<Jn;j++){
      temp1 = arma::zeros<arma::mat>(p+1, p+1);
      temp_dvalue = 0;
      for(l=0;l<n;l++){
        Kn = KerexpgamZ_list(i*n + l).n_cols;
        for(k=0;k<Kn;k++){
          temp_vec1 = Xmat_list(i*n+l).tube(j,k);
          temp_vec2 = KerXexpgamZ_list(i*n + l).tube(j,k);
          temp1 = temp1 + temp_vec2 * temp_vec1.t();
        }
        temp_dvalue = temp_dvalue + arma::accu(KerexpgamZ_list(i*n + l).row(j));
      }
      XXtbar_list(i).slice(j) = temp1 / temp_dvalue;
    }
    //XXtbar_list(i).transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
  }

  // TODO fix outerker_Y_C, use outer input, change the type of temp1 and temp0, they are vec
  // Calculate \bar{Y}
  // arma::field<arma::vec> Ybar_list(n);
  // arma::vec temp_response_i,temp_response_l ;
  // arma::mat temp_kermat_T;
  // for(i=0;i<n;i++){
  //   temp_response_i = vec(response[i].begin(),response[i].size(),false);
  //   temp1 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, 1);
  //   temp0 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, 1);
  //
  //
  //   for(l=0;l<n;l++){
  //     temp_response_l = vec(response[l].begin(),response[l].size(),false);
  //     temp_kermat_T = outerker_Y_C(temp_response_i,temp_response_l,pow(n,-0.8));
  //     temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
  //
  //     //temp = arma::sum(KerexpgamZ_list(i*n + l),1);
  //     temp = arma::sum(temp_kermat,1);
  //     temp_vec1 = temp_kermat_T * temp_response_l;
  //     //temp_vec1.transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
  //     temp1 = temp1 + temp % temp_vec1 ;
  //     temp0 = temp0 + temp % arma::sum(temp_kermat_T,1);
  //
  //   }
  //   Ybar_list(i) = temp1 /temp0;
  //   for(j=0;j<Ybar_list(i).size();j++){
  //     if(std::isnan(Ybar_list(i).at(j))){
  //       Ybar_list(i).at(j) = temp_response_i.at(j);
  //     }
  //   }
  //
    //Rcpp::Rcout << temp0 <<endl;
    //Ybar_list(i).transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
  //}



  // Calculate theta
  //Rcpp::Rcout <<  "3step" <<endl;


  temp1 = arma::zeros<arma::vec>(p+1);// numerator
  temp0 = arma::zeros<arma::mat>(p+1,p+1); // denominator
  arma::vec temp_dlambda,temp_response;
  // for(i=0;i<n;i++){
  //
  //   temp_dlambda = vec(dlambda[i].begin(),dlambda[i].size(),false);
  //   for(l=0;l<n;l++){
  //
  //     Jn = KerexpgamZ_list(i*n + l).n_rows;
  //     Kn = KerexpgamZ_list(i*n + l).n_cols;
  //
  //     for(j=0;j<Jn;j++){
  //       for(k=0;k<Kn;k++){
  //         temp_vec1 = KerXexpgamZ_list(i*n+l).tube(j,k);
  //         temp_vec2 = Xmat_list(i*n+l).tube(j,k);
  //         temp_vec2 = temp_vec2-Xbar_list(i).row(j).t();
  //         temp0 = temp0+temp_dlambda(j)*temp_vec1*temp_vec2.t();
  //       }
  //     }
  //   }
  // }
  //
  // temp_test = temp0;
  //
  // Rcpp::Rcout << temp0 << endl;

  temp0 = arma::zeros<arma::mat>(p+1,p+1); // denominator
  for(i=0;i<n;i++){
      Jn = KerexpgamZ_list(i*n + i).n_rows;
      Kn = KerexpgamZ_list(i*n + i).n_cols;
      for(j=0;j<Jn;j++){
        for(k=0;k<Kn;k++){
          temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
          //temp = Xbar_list(i).row(j);
          temp0 = temp0 + temp_kermat(j,k)*(XXtbar_list(i).slice(j) - Xbar_list(i).row(j).t() * Xbar_list(i).row(j));
        }
      }
  }

  //Rcpp::Rcout << temp0 << endl;
  //Rcpp::Rcout <<  "4step" <<endl;
  for(i=0;i<n;i++){
    Jn = KerexpgamZ_list(i*n + i).n_rows;
    Kn = KerexpgamZ_list(i*n + i).n_cols;
    temp_response = vec(response[i].begin(),response[i].size(),false);
    temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
    temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
    for(j=0;j<Jn;j++){
      for(k=0;k<Kn;k++){
        temp_vec2 = Xmat_list(i*n+i).tube(j,k);
        //temp1 = temp1+(temp_response(j)-Ybar_list(i).at(j))*(temp_vec2-Xbar_list(i).row(j).t()) * temp_kermat.at(j,k) ;
        //temp1 = temp1+(temp_response(j)*temp_vec2-Ybar_list(i).at(j)*Xbar_list(i).row(j).t()) * temp_kermat.at(j,k) ;
        temp1 = temp1+temp_response(j)*(temp_vec2-Xbar_list(i).row(j).t()) * temp_kermat.at(j,k);
      }
    }
    //temp1.t().print();
    //if(temp1.has_nan()) {
    //  Rcpp::Rcout << i <<  endl;
    //}
  }
  arma::vec  thetaest= inv_sympd(temp0) * temp1;
  //Rcpp::Rcout << inv_sympd(temp_test) * temp1 << endl;
  // return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
  //                           Rcpp::Named("expgamZ_list") =  expgamZ_list,
  //                           Rcpp::Named("KerexpgamZ_list") = KerexpgamZ_list,
  //                           Rcpp::Named("Xmat_list") = Xmat_list,
  //                           Rcpp::Named("KerXexpgamZ_list") = KerXexpgamZ_list,
  //                           Rcpp::Named("Xbar_list") = Xbar_list,
  //                           Rcpp::Named("temp0") = temp0,
  //                           Rcpp::Named("temp1") = temp1
  //                         );
  return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
                            Rcpp::Named("Xbar_list") = Xbar_list,
                            Rcpp::Named("temp0") = temp0,
                            Rcpp::Named("temp1") = temp1

  );
}



// Rcpp::List longest_test_c(const arma::rowvec& gamma,
//                      Rcpp::ListOf<NumericMatrix>& kerMat,
//                      Rcpp::ListOf<NumericVector>& meas_times,
//                      Rcpp::ListOf<NumericMatrix>& covariates,
//                      Rcpp::ListOf<NumericVector>& response,
//                      Rcpp::ListOf<NumericVector>& dlambda,
//                      const arma::vec& censor, const unsigned int& n,
//                      const unsigned int& p) {
//   unsigned int i, l=0;
//   arma::mat temp_kermat,temp_test;
//   arma::mat temp1, temp0, s0sumres, s1sumres;
//   arma::cube temp_cube;
//   arma::mat tempres;
//   double temp_dvalue;
//   arma::dmat temp_cov_i;
//   arma::rowvec expgammaZ;
//   arma::field<arma::rowvec> expgamZ_list(n);
//   arma::field<arma::mat> Xbar_list(n);
//   arma::field<arma::mat> KerexpgamZ_list(n*n);
//   arma::field<arma::cube> Xmat_list(n*n),KerXexpgamZ_list(n*n), XXtbar_list(n);
//   arma::vec temp_meas_time_i,temp_meas_time_l,temp_countprocess_l,temp_vec1,temp_vec2;
//   arma::vec censorind;
//
//   unsigned int j,k,Jn,Kn;
//   //Rcpp::Rcout <<  "start" <<endl;
//   //List res(2);
//   // i is the iterator for R_i
//   for (i = 0; i < n; i++) {
//     temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
//     temp_cov_i = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
//     // Obtain exp(gamma * Z_i(R_{ik}))
//     // Store them in expgamZ_list
//
//     expgamZ_list(i) = conv_to<arma::rowvec>::from(exp(gamma * temp_cov_i));
//   }
//
//   //Rcpp::Rcout << expgamZ_list <<endl;
//   //Rcpp::Rcout << KerexpgamZ_list <<endl;
//   //Rcpp::Rcout <<  "2step" <<endl;
//
//   // Calculate \bar{X}
//   arma::vec tempres_vec;
//   arma::mat temp_cov;
//   arma::vec temp = arma::vec(p+1);
//   for(i=0;i<n;i++){
//     //temp1 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, p+1);
//
//     temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
//     temp_kermat = mat(kerMat[i * n + 0].begin(),kerMat[i * n + 0].nrow(),kerMat[i * n + 0].ncol(),false);
//     Jn = temp_kermat.n_rows;
//     Xbar_list(i) = arma::mat(p+1,Jn);
//     for(j=0;j<Jn;j++){
//       temp_dvalue = 0;
//       tempres_vec = arma::zeros<arma::vec>(p+1);
//       for(l=0;l<n;l++){
//         temp_meas_time_l = vec(meas_times[l].begin(),meas_times[l].size(),false);
//         temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
//         temp_countprocess_l = countprofun_C(temp_meas_time_l, temp_meas_time_i);
//         temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
//         Kn = temp_kermat.n_cols;
//         for(k=0;k<Kn;k++){
//           temp.head_rows(p) = temp_cov.col(k);
//           temp(p) = temp_countprocess_l(j);
//           temp_dvalue = temp_dvalue + temp_kermat(j,k) * (censor[l] > temp_meas_time_i(j)) * expgamZ_list(l).at(k);
//           tempres_vec = tempres_vec + temp_kermat(j,k) * temp * (censor[l] > temp_meas_time_i(j)) * expgamZ_list(l).at(k);
//           //Rcpp::Rcout <<i << j << l << k << endl;
//         }
//       }
//       Xbar_list(i).col(j) = tempres_vec/temp_dvalue ;
//     }
//     Xbar_list(i).transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
//   }
//
//   // TODO futher explain
//   // Calculate \bar{XXt}
//   //Rcpp::Rcout <<  "2.5step" <<endl;
//   arma::mat tempres_mat;
//   for(i=0;i<n;i++){
//     //temp1 = arma::zeros<arma::mat>(KerexpgamZ_list(i*n+0).n_rows, p+1);
//
//     temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
//     temp_kermat = mat(kerMat[i * n + 0].begin(),kerMat[i * n + 0].nrow(),kerMat[i * n + 0].ncol(),false);
//     Jn = temp_kermat.n_rows;
//     XXtbar_list(i) = arma::cube(p+1,p+1,Jn);
//     for(j=0;j<Jn;j++){
//       temp_dvalue = 0;
//       tempres_mat = arma::zeros<arma::mat>(p+1,p+1);
//       for(l=0;l<n;l++){
//         temp_meas_time_l = vec(meas_times[l].begin(),meas_times[l].size(),false);
//         temp_cov = mat(covariates[l].begin(),covariates[l].nrow(),covariates[l].ncol(),false);
//         temp_countprocess_l = countprofun_C(temp_meas_time_l, temp_meas_time_i);
//         temp_kermat = mat(kerMat[i * n + l].begin(),kerMat[i * n + l].nrow(),kerMat[i * n + l].ncol(),false);
//         Kn = temp_kermat.n_cols;
//         for(k=0;k<Kn;k++){
//           temp.head_rows(p) = temp_cov.col(k);
//           temp(p) = temp_countprocess_l(j);
//           temp_dvalue = temp_dvalue + temp_kermat(j,k) * (censor[l] > temp_meas_time_i(j)) * expgamZ_list(l).at(k);
//           tempres_mat = tempres_mat + temp_kermat(j,k) * temp * temp .t() * (censor[l] > temp_meas_time_i(j)) * expgamZ_list(l).at(k);
//
//         }
//       }
//       XXtbar_list(i).slice(j) = tempres_mat / temp_dvalue;
//     }
//     XXtbar_list(i).transform( [](double val) { return (std::isnan(val) ? 0 : val); } );
//     //Rcpp::Rcout << XXtbar_list(i) << endl;
//   }
//
//
//   // Calculate theta
//   //Rcpp::Rcout <<  "3step" <<endl;
//
//
//   temp1 = arma::zeros<arma::vec>(p+1);// numerator
//   temp0 = arma::zeros<arma::mat>(p+1,p+1); // denominator
//   arma::vec temp_dlambda,temp_response;
//   // for(i=0;i<n;i++){
//   //
//   //   temp_dlambda = vec(dlambda[i].begin(),dlambda[i].size(),false);
//   //   for(l=0;l<n;l++){
//   //
//   //     Jn = KerexpgamZ_list(i*n + l).n_rows;
//   //     Kn = KerexpgamZ_list(i*n + l).n_cols;
//   //
//   //     for(j=0;j<Jn;j++){
//   //       for(k=0;k<Kn;k++){
//   //         temp_vec1 = KerXexpgamZ_list(i*n+l).tube(j,k);
//   //         temp_vec2 = Xmat_list(i*n+l).tube(j,k);
//   //         temp_vec2 = temp_vec2-Xbar_list(i).row(j).t();
//   //         temp0 = temp0+temp_dlambda(j)*temp_vec1*temp_vec2.t();
//   //       }
//   //     }
//   //   }
//   // }
//   //
//   // temp_test = temp0;
//   //
//   // Rcpp::Rcout << temp0 << endl;
//   //Rcpp::Rcout << "Step 4" << endl;
//   //temp0 = arma::zeros<arma::mat>(p+1,p+1); // denominator
//   for(i=0;i<n;i++){
//     temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n +i].nrow(),kerMat[i * n + i].ncol(),false);
//     Jn = temp_kermat.n_rows;
//     Kn = temp_kermat.n_cols;
//     //Rcpp::Rcout <<Jn <<" " << Kn << endl;
//     for(j=0;j<Jn;j++){
//       for(k=0;k<Kn;k++){
//         //temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
//         tempres_vec = Xbar_list(i).col(j);
//         temp0 = temp0 + temp_kermat(j,k)*(XXtbar_list(i).slice(j) - tempres_vec * tempres_vec.t());
//         //Rcpp::Rcout << (XXtbar_list(i).slice(j) - tempres_vec * tempres_vec.t()) << endl;
//       }
//     }
//
//   }
//
//   //Rcpp::Rcout << temp0 << endl;
//   //Rcpp::Rcout <<  "4step" <<endl;
//   for(i=0;i<n;i++){
//     temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n +i].nrow(),kerMat[i * n + i].ncol(),false);
//     Jn = temp_kermat.n_rows;
//     Kn = temp_kermat.n_cols;
//     temp_cov = mat(covariates[i].begin(),covariates[i].nrow(),covariates[i].ncol(),false);
//     temp_response = vec(response[i].begin(),response[i].size(),false);
//     //temp_kermat = mat(kerMat[i * n + i].begin(),kerMat[i * n + i].nrow(),kerMat[i * n + i].ncol(),false);
//     temp_meas_time_i = vec(meas_times[i].begin(),meas_times[i].size(),false);
//
//     for(j=0;j<Jn;j++){
//       for(k=0;k<Kn;k++){
//         temp.head_rows(p) = temp_cov.col(k);
//         temp(p) = j;
//         temp1 = temp1+temp_response(j)*(temp-Xbar_list(i).col(j)) * temp_kermat.at(j,k) ;
//       }
//     }
//
//   }
//   arma::vec  thetaest= inv_sympd(temp0) * temp1;
//   //Rcpp::Rcout << inv_sympd(temp_test) * temp1 << endl;
//   // return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
//   //                           Rcpp::Named("expgamZ_list") =  expgamZ_list,
//   //                           Rcpp::Named("KerexpgamZ_list") = KerexpgamZ_list,
//   //                           Rcpp::Named("Xmat_list") = Xmat_list,
//   //                           Rcpp::Named("KerXexpgamZ_list") = KerXexpgamZ_list,
//   //                           Rcpp::Named("Xbar_list") = Xbar_list,
//   //                           Rcpp::Named("temp0") = temp0,
//   //                           Rcpp::Named("temp1") = temp1
//   //                         );
//   return Rcpp::List::create(Rcpp::Named("thetaest") = thetaest,
//                             Rcpp::Named("temp0") = temp0,
//                             Rcpp::Named("temp1") = temp1);
// }






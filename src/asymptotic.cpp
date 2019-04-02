#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]

#include <iostream>

#include "longest.hpp"

using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
Rcpp::List Xbar_c(const arma::rowvec &gamma,
                  const arma::rowvec &theta,
                  Rcpp::ListOf<NumericMatrix> &kerMat,
                  Rcpp::ListOf<NumericVector> &meas_times,
                  Rcpp::ListOf<NumericMatrix> &covariates,
                  const arma::vec &censor,
                  const unsigned int &n,
                  const unsigned int &p)
{
  unsigned int i, l = 0;
  arma::mat temp_kermat, temp_cov;
  arma::mat temp_Xbar_list;
  arma::cube temp_XXtbar_list, temp_XXtZbar_list;
  arma::vec temp_S0_list;

  arma::vec Xbar_temp, X_temp = vec(p + 1);
  arma::mat XXtbar_temp = mat(p + 1, p + 1);
  arma::mat XXtZbar_temp = mat(p + 1, p);

  arma::rowvec expgammaZ;
  //arma::field<arma::rowvec> expgamZ_list(n);

  // Store list of results to return
  Rcpp::List Xbar_list(n);
  Rcpp::List XXtbar_list(n);
  Rcpp::List XXtZbar_list(n);
  Rcpp::List S0_list(n);
  //arma::field<arma::mat> KerexpgamZ_list(n * n);

  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;
  arma::vec censorind, temp_response;
  double weight;
  arma::uvec temp_uvec;

  arma::uvec temp_counting;
  double S0_temp;

  unsigned int j, k, Jn, Kn;

  for (i = 0; i < n; i++)
  {
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    Jn = temp_meas_time_i.size();
    temp_Xbar_list = mat(Jn, p + 1, arma::fill::zeros);
    temp_XXtbar_list = cube(p + 1, p + 1, Jn, arma::fill::zeros);
    temp_XXtZbar_list = cube(p + 1, p, Jn, arma::fill::zeros);
    temp_S0_list = vec(Jn, arma::fill::zeros);

    for (j = 0; j < Jn; j++)
    {
      Xbar_temp = vec(p + 1, arma::fill::zeros);
      S0_temp = 0;
      XXtbar_temp = mat(p + 1, p + 1, arma::fill::zeros);
      XXtZbar_temp = mat(p + 1, p, arma::fill::zeros);

      for (l = 0; l < n; l++)
      {
        if (censor[l] > temp_meas_time_i[j])
        {
          temp_meas_time_l =
            vec(meas_times[l].begin(), meas_times[l].size(), false);
          temp_cov = mat(covariates[l].begin(), covariates[l].nrow(),
                         covariates[l].ncol(), false);
          temp_kermat = mat(kerMat[i * n + l].begin(), kerMat[i * n + l].nrow(),
                            kerMat[i * n + l].ncol(), false);
          Kn = temp_cov.n_cols;

          for (k = 0; k < Kn; k++)
          {
            if (abs(temp_kermat(j, k)) > 1e-10)
            {
              X_temp.head(p) = temp_cov.col(k);
              temp_uvec = find(temp_meas_time_l >= temp_meas_time_i(j), 1);

              X_temp(p) = temp_uvec.n_elem == 0 ? temp_meas_time_l.n_elem : temp_uvec(0);
              // cout << temp_kermat(j, k) * exp(gamma * temp_cov.col(k)) <<endl;
              weight = as_scalar(temp_kermat(j, k) * exp(gamma * temp_cov.col(k)));

              Xbar_temp = Xbar_temp + X_temp * weight;
              XXtbar_temp = XXtbar_temp + X_temp * X_temp.t() * weight;
              XXtZbar_temp = XXtZbar_temp + X_temp * X_temp.t() * theta.t() * X_temp.head(p).t() * weight;
              S0_temp = S0_temp + weight;
            }
          }
        }
      }
      temp_S0_list(j) = S0_temp;

      if (S0_temp != 0)
      {
        temp_Xbar_list.row(j) = trans(Xbar_temp / S0_temp);

        temp_XXtbar_list.slice(j) = XXtbar_temp / S0_temp;

        temp_XXtZbar_list.slice(j) = XXtZbar_temp / S0_temp;
      }
    }
    Xbar_list(i) = temp_Xbar_list;
    S0_list(i) = temp_S0_list;
    XXtbar_list(i) = temp_XXtbar_list;
    XXtZbar_list(i) = temp_XXtZbar_list;
  }
  return Rcpp::List::create(Rcpp::Named("Xbar") = Xbar_list,
                            Rcpp::Named("XXtbar") = XXtbar_list,
                            Rcpp::Named("XXtZbar") = XXtZbar_list,
                            Rcpp::Named("S0") = S0_list);
}

// // [[Rcpp::export]]
// Rcpp::List Rs_c(const arma::rowvec & gamma,
//                 const arma::rowvec & theta,
//                 Rcpp::ListOf < NumericMatrix > & kerMat,
//                 Rcpp::ListOf < NumericVector > & meas_times,
//                 Rcpp::ListOf < NumericMatrix > & covariates,
//                 Rcpp::ListOf < NumericMatrix > & Xbar,
//                 Rcpp::ListOf < NumericVector > & S0,
//                 const arma::vec & censor,
//                 const unsigned int & n,
//                 const unsigned int & p) {
//   unsigned int i, l = 0;
//   arma::mat temp_kermat, temp_cov, temp_Rs_list, Xbar_temp;
//   arma::vec X_temp, Rs_temp = vec(p + 1);
//   arma::rowvec expgammaZ;
//   arma::field < arma::rowvec > expgamZ_list(n);
//   Rcpp::List Rs_list(n);
//   arma::field < arma::mat > KerexpgamZ_list(n * n);
//
//   arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;
//   arma::vec censorind, temp_S0;
//   double weight;
//   arma::uvec temp_uvec;
//
//   arma::uvec temp_counting;
//
//   unsigned int j, k, Jn, Kn;
//   X_temp = vec(p + 1, arma::fill::zeros);
//   for (i = 0; i < n; i++) {
//     //cout << "i: " << i << endl;
//     temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
//     Jn = temp_meas_time_i.size();
//     temp_Rs_list = mat(Jn, p + 1, arma::fill::zeros);
//     temp_S0 = vec(S0[i].begin(), S0[i].size(), false);
//
//     for (j = 0; j < Jn; j++) {
//       //cout << "j: " << j << endl;
//       if (temp_S0(j) != 0) {
//         Rs_temp = vec(p + 1, arma::fill::zeros);
//
//         Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
//                         Xbar[i].ncol(), false);
//
//         for (l = 0; l < n; l++) {
//           //cout << "  l: " << l << endl;
//           if (censor[l] > temp_meas_time_i[j]) {
//             temp_meas_time_l =
//               vec(meas_times[l].begin(), meas_times[l].size(), false);
//             temp_cov = mat(covariates[l].begin(), covariates[l].nrow(),
//                            covariates[l].ncol(), false);
//             temp_kermat = mat(kerMat[i * n + l].begin(), kerMat[i * n + l].nrow(),
//                               kerMat[i * n + l].ncol(), false);
//             Kn = temp_cov.n_cols;
//
//             for (k = 0; k < Kn; k++) {
//               //cout << "    k: " << k << "/" << Kn << endl;
//               if (abs(temp_kermat(j, k)) > 1e-10) {
//
//                 X_temp.head(p) = temp_cov.col(k);
//                 temp_uvec = find(temp_meas_time_l >= temp_meas_time_i(j), 1);
//
//                 X_temp(p) = temp_uvec.n_elem == 0 ? temp_meas_time_l.n_elem : temp_uvec(0);
//
//                 weight = as_scalar(temp_kermat(j, k) * exp(gamma * temp_cov.col(k)));
//
//                 Rs_temp = Rs_temp + (X_temp - trans(Xbar_temp.row(j))) * as_scalar(theta * X_temp) * weight;
//
//               }
//             }
//           }
//         }
//         //cout << "Rs_temp: " << endl;
//         //Rs_temp.print();
//         //cout << "S0_temp: " << temp_S0.t() << endl;
//
//         temp_Rs_list.row(j) = trans(Rs_temp / temp_S0(j));
//       }
//     }
//     Rs_list(i) = temp_Rs_list;
//
//   }
//   return Rs_list;
// }

// [[Rcpp::export]]
Rcpp::List H_A_c(const arma::rowvec &gamma,
                 const arma::rowvec &theta,
                 Rcpp::ListOf<NumericMatrix> &kerMat,
                 Rcpp::ListOf<NumericVector> &meas_times,
                 Rcpp::ListOf<NumericMatrix> &covariates,
                 Rcpp::ListOf<NumericVector> &response,
                 Rcpp::ListOf<NumericMatrix> &Xbar,
                 Rcpp::ListOf<NumericVector> &XXtbar,
                 Rcpp::ListOf<NumericVector> &XXtZbar,
                 Rcpp::ListOf<NumericVector> &S0,
                 const arma::vec &censor,
                 const unsigned int &n,
                 const unsigned int &p)
{
  unsigned int i = 0;
  //Rcpp::Rcout << "0" << endl;
  arma::mat temp_kermat, temp_cov, temp_Rs_list;
  arma::vec X_temp, temp_X_diff = vec(p + 1);
  arma::mat A_temp = mat(p, p, arma::fill::zeros);
  arma::mat H_temp = mat(p + 1, p, arma::fill::zeros);

  arma::rowvec expgammaZ;
  arma::field<arma::rowvec> expgamZ_list(n);

  arma::field<arma::mat> KerexpgamZ_list(n * n);

  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l, temp_response;
  arma::vec censorind, temp_S0;
  arma::mat Xbar_temp;
  arma::cube XXtbar_temp, XXtZbar_temp;
  arma::mat Xbar_dev_temp;
  arma::vec S0_temp;

  unsigned int j, k, Jn, Kn;
  //X_temp = vec(p + 1, arma::fill::zeros);
  for (i = 0; i < n; i++)
  {

    //Rcpp::Rcout << "1" << endl;

    temp_response = vec(response[i].begin(), response[i].size(), false);

    Jn = temp_response.size();
    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);
    temp_cov = mat(covariates[i].begin(), covariates[i].nrow(),
                   covariates[i].ncol(), false);
    temp_kermat = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                      kerMat[i * n + i].ncol(), false);

    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);

    //Rcpp::Rcout << "2" << endl;

    XXtbar_temp = arma::cube(XXtbar[i].begin(), p + 1, p + 1, Jn, false);

    //XXtbar_temp.print();
    //Rcpp::Rcout << "3" << endl;

    XXtZbar_temp = arma::cube(XXtZbar[i].begin(), p + 1, p, Jn, false);

    S0_temp = vec(S0[i].begin(), S0[i].size(), false);


    for (j = 0; j < Jn; j++)
    {

      Kn = temp_cov.n_cols;

      //cout << Xbar_temp.row(j).head(p).t()<<endl;

      //cout << XXtbar_temp.slice(j)<<endl;

      for (k = 0; k < Kn; k++)
      {
        if (abs(temp_kermat(j, k)) > 1e-10)
        {

          // Calculate H
          //Rcpp::Rcout << "H" << endl;
          Xbar_dev_temp = XXtbar_temp.slice(j).submat(0, 0, p, p - 1) -
            Xbar_temp.row(j).t() * Xbar_temp.row(j).head(p);

          H_temp = H_temp - temp_kermat(j, k) * temp_response(j) * Xbar_dev_temp;

          H_temp = H_temp - temp_kermat(j, k) * (XXtZbar_temp.slice(j) -
            XXtbar_temp.slice(j) * theta.t() * Xbar_temp.row(j).head(p));

          //H_temp.t().print();

          H_temp = H_temp + temp_kermat(j, k)*(as_scalar(Xbar_temp.row(j) * theta.t()) * Xbar_dev_temp+
            Xbar_temp.row(j).t() * theta * Xbar_dev_temp);

          //H_temp.t().print();

          A_temp = A_temp + temp_kermat(j, k) * (XXtbar_temp.slice(j).submat(0, 0, p - 1, p - 1) -
            Xbar_temp.row(j).head(p).t() * Xbar_temp.row(j).head(p));


        }
      }
    }

  }
  return Rcpp::List::create(Rcpp::Named("Hmat") = H_temp,
                            Rcpp::Named("Amat") = A_temp);
}

// [[Rcpp::export]]
Rcpp::List long_asy_c(const arma::rowvec &gamma,
                      const arma::rowvec &theta,
                      Rcpp::ListOf<NumericMatrix> &kerMat,
                      Rcpp::ListOf<NumericVector> &meas_times,
                      Rcpp::ListOf<NumericMatrix> &covariates,
                      Rcpp::ListOf<NumericMatrix> &Xbar,
                      Rcpp::ListOf<NumericVector> &XXtbar,
                      Rcpp::ListOf<NumericVector> &response,
                      const arma::mat &Hmat,
                      const arma::mat &Amat,
                      const arma::vec &censor,
                      const unsigned int &n,
                      const unsigned int &p)
{

  //cout << "HHH" <<endl;
  unsigned int i = 0;
  arma::mat temp_kermat, temp_cov, Xbar_temp, temp_response;
  arma::vec X_temp, temp_part1, temp_part2, temp_part3 = vec(p + 1);
  arma::mat temp_X_diff = mat(p + 1, p + 1);
  arma::cube XXtbar_temp;
  //arma::mat D_temp = mat(p,p,arma::fill::zeros);
  //arma::mat P_temp = mat(p+1,p,arma::fill::zeros);
  arma::rowvec expgammaZ;
  //arma::field<arma::rowvec> expgamZ_list(n);
  //Rcpp::List Rs_list(n);
  //arma::field<arma::mat> KerexpgamZ_list(n * n);

  arma::mat V_temp = mat(p, p, arma::fill::zeros);

  arma::mat res_temp_ind = vec(p + 1);
  arma::mat Bmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat Sigmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l,gmu_temp, temp_lambda ;

  arma::vec censorind, temp_S0;

  arma::uvec temp_uvec;

  arma::uvec temp_counting;

  arma::mat HAi = Hmat * inv(Amat);

  unsigned int j, k, Jn, Kn;
  X_temp = vec(p + 1, arma::fill::zeros);

  temp_part3 = vec(p, arma::fill::zeros);

  arma::vec zdiff = vec(n,arma::fill::zeros);

  double weight;
  for (i = 0; i < n; i++)
  {
    //cout << "i: " << i << endl;
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    temp_response = vec(response[i].begin(), response[i].size(), false);
    // temp_lambda = vec(asylambda[i].begin(), lambda[i].size(), false);
    // gmu_temp = vec(gmu[i].begin(), gmu[i].size(), false);
    //temp_dlambda = vec(dlambda[i].begin(), dlambda[i].size(), false);
    //temp_mu0 = vec(mu0[i].begin(), mu0[i].size(), false);

    //temp_Rs = vec(Rs[i].begin(), Rs[i].size(), false);

    temp_cov = mat(covariates[i].begin(), covariates[i].nrow(),
                   covariates[i].ncol(), false);
    temp_kermat = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                      kerMat[i * n + i].ncol(), false);

    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);

    Jn = temp_meas_time_i.size();
    Kn = temp_cov.n_cols;

    XXtbar_temp = arma::cube(XXtbar[i].begin(), p + 1, p + 1, Jn, false);

    temp_part1 = vec(p + 1, arma::fill::zeros);
    temp_part2 = vec(p, arma::fill::zeros);






    for (j = 0; j < Jn; j++)
    {

      //cout << "j: " << j << endl;

      for (k = 0; k < Kn; k++)
      {
        //cout << " k: " << k << "/" << Kn << endl;
        //cout << "Kernal" << temp_kermat(j, k) << endl;

        if (abs(temp_kermat(j, k)) > 1e-10)
        {

          X_temp.head(p) = temp_cov.col(k);
          //temp_uvec = find(temp_meas_time_l >= temp_meas_time_i(j),1);

          X_temp(p) = j;

          //cout << "X_temp: " << X_temp.t() <<endl;

          //cout << "XXtbar_tem: " << XXtbar_temp.slice(j) <<endl;

          temp_X_diff = XXtbar_temp.slice(j) - Xbar_temp.row(j).t() * Xbar_temp.row(j);

          Bmat = Bmat - temp_kermat(j, k) * temp_X_diff;

          //temp_X_diff = X_temp - trans(Xbar_temp.row(j));
          ////================= PART 1 ===================////
          //cout << "PART1";
          //weight = 1-as_scalar(exp(gamma * X_temp.head(p)))*temp_lambda.at(j);
          weight = 1;


          temp_part1 = temp_part1 + temp_kermat(j, k) * (temp_response(j) * (X_temp - Xbar_temp.row(j).t()) -
            temp_X_diff * theta.t())*weight;

          //temp_part1 = temp_part1 + temp_kermat(j, k) *(temp_response(j)-(temp_response(j)+
          //  as_scalar(theta * (X_temp - Xbar_temp.row(j).t())))*as_scalar(exp(gamma * X_temp.head(p)))*temp_lambda.at(j));

          temp_part2 = temp_part2 + temp_kermat(j, k) * (X_temp.head(p) - Xbar_temp.row(j).head(p).t())*weight;

          // Rcpp::Rcout << "i: " << setw(2) << i << " j: " << setw(2) << j <<
          //   " k: " << setw(2) << k << " ker: " << temp_kermat(j, k);
          // Rcpp::Rcout << " part2: " << weight<<  endl;




        }
      }
      //cout << "end" << endl;
    }



    //cout << temp_part1 << endl;
    //cout << temp_part2 << endl;

    res_temp_ind = temp_part1 + HAi * temp_part2;
    V_temp = V_temp + temp_part2 * temp_part2.t();
    //zdiff(i) = as_scalar(temp_part2);
    //temp_part3 += temp_part2 * temp_part2;
    //res_temp_ind  = temp_part1;
    //cout << "Z: ";
    //cout << temp_part2.t() << endl;
    Sigmat = Sigmat + res_temp_ind * res_temp_ind.t();
  }
  //cout <<"============:" << temp_part3 << endl;
  Sigmat = Sigmat;
  Bmat = Bmat;
  V_temp = V_temp;
  return Rcpp::List::create(Rcpp::Named("Sigmat") = Sigmat,
                            Rcpp::Named("Bmat") = Bmat,
                            Rcpp::Named("Vmat") = V_temp);
}

// [[Rcpp::export]]
arma::mat EE_c(const arma::rowvec &gamma,
               const arma::rowvec &theta,
               Rcpp::ListOf<NumericMatrix> &kerMat,
               Rcpp::ListOf<NumericVector> &meas_times,
               Rcpp::ListOf<NumericMatrix> &covariates,
               Rcpp::ListOf<NumericMatrix> &Xbar,
               Rcpp::ListOf<NumericVector> &XXtbar,
               Rcpp::ListOf<NumericVector> &response,
               const arma::mat &Hmat,
               const arma::mat &Amat,
               const arma::vec &censor,
               const unsigned int &n,
               const unsigned int &p)
{

  unsigned int i = 0;
  arma::mat temp_kermat, temp_cov, Xbar_temp, temp_response;
  arma::vec X_temp, temp_part1, temp_part2, temp_part3 = vec(p + 1);
  arma::mat temp_X_diff = mat(p + 1, p + 1);
  arma::cube XXtbar_temp;

  arma::rowvec expgammaZ;

  arma::mat res_temp_ind = vec(p + 1);
  arma::mat Bmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat Sigmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat res = mat(p + 1, 1, arma::fill::zeros);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;

  arma::vec censorind, temp_S0;
  arma::uvec temp_uvec;
  arma::uvec temp_counting;
  arma::mat HAi = Hmat * inv(Amat);

  unsigned int j, k, Jn, Kn;
  X_temp = vec(p + 1, arma::fill::zeros);
  for (i = 0; i < n; i++)
  {
    //cout << "i: " << i << endl;
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    temp_response = vec(response[i].begin(), response[i].size(), false);
    temp_cov = mat(covariates[i].begin(), covariates[i].nrow(),
                   covariates[i].ncol(), false);
    temp_kermat = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                      kerMat[i * n + i].ncol(), false);

    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);

    Jn = temp_meas_time_i.size();
    Kn = temp_cov.n_cols;

    XXtbar_temp = arma::cube(XXtbar[i].begin(), p + 1, p + 1, Jn, false);

    temp_part1 = vec(p + 1, arma::fill::zeros);
    temp_part2 = vec(p, arma::fill::zeros);

    for (j = 0; j < Jn; j++)
    {
      for (k = 0; k < Kn; k++)
      {
        if (abs(temp_kermat(j, k)) > 0)
        {
          X_temp.head(p) = temp_cov.col(k);
          X_temp(p) = j;
          temp_X_diff = XXtbar_temp.slice(j) - Xbar_temp.row(j).t() * Xbar_temp.row(j);

          Bmat = Bmat - temp_kermat(j, k) * temp_X_diff;
          temp_part1 = temp_part1 + temp_kermat(j, k)*(temp_response(j)*(X_temp - Xbar_temp.row(j).t())-
             temp_X_diff*theta.t());
          //temp_part1 = temp_part1 + temp_kermat(j, k) * (Xbar_temp.row(j).t() * Xbar_temp.row(j)*theta.t());
        }
      }
      //cout << "end" << endl;
    }

    //cout << temp_part1 << endl;
    //cout << temp_part2 << endl;
    res = res + temp_part1;
    //res_temp_ind  = temp_part1;
  }
  return res;
}

// [[Rcpp::export]]
arma::mat EEZ_c(const arma::rowvec &gamma,
                const arma::rowvec &theta,
                Rcpp::ListOf<NumericMatrix> &kerMat,
                Rcpp::ListOf<NumericVector> &meas_times,
                Rcpp::ListOf<NumericMatrix> &covariates,
                Rcpp::ListOf<NumericMatrix> &Xbar,
                Rcpp::ListOf<NumericVector> &XXtbar,
                Rcpp::ListOf<NumericVector> &response,
                const arma::mat &Hmat,
                const arma::mat &Amat,
                const arma::vec &censor,
                const unsigned int &n,
                const unsigned int &p)
{

  unsigned int i = 0;
  arma::mat temp_kermat, temp_cov, Xbar_temp, temp_response;
  arma::vec X_temp, temp_part1, temp_part2, temp_part3 = vec(p + 1);
  arma::mat temp_X_diff = mat(p + 1, p + 1);
  arma::cube XXtbar_temp;

  arma::rowvec expgammaZ;

  arma::mat res_temp_ind = vec(p + 1);
  arma::mat Bmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat Sigmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat res = mat(p, 1, arma::fill::zeros);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;

  arma::vec censorind, temp_S0;

  arma::uvec temp_uvec;

  arma::uvec temp_counting;

  arma::mat HAi = Hmat * inv(Amat);

  unsigned int j, k, Jn, Kn;
  X_temp = vec(p + 1, arma::fill::zeros);
  for (i = 0; i < n; i++)
  {
    //cout << "i: " << i << endl;
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);
    temp_response = vec(response[i].begin(), response[i].size(), false);

    temp_cov = mat(covariates[i].begin(), covariates[i].nrow(),
                   covariates[i].ncol(), false);
    temp_kermat = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                      kerMat[i * n + i].ncol(), false);

    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);

    Jn = temp_meas_time_i.size();
    Kn = temp_cov.n_cols;

    XXtbar_temp = arma::cube(XXtbar[i].begin(), p + 1, p + 1, Jn, false);

    temp_part1 = vec(p + 1, arma::fill::zeros);
    temp_part2 = vec(p, arma::fill::zeros);

    for (j = 0; j < Jn; j++)
    {
      for (k = 0; k < Kn; k++)
      {

        if (abs(temp_kermat(j, k)) > 0)
        {

          X_temp.head(p) = temp_cov.col(k);
          X_temp(p) = j;

          temp_X_diff = XXtbar_temp.slice(j) - Xbar_temp.row(j).t() * Xbar_temp.row(j);

          //Bmat = Bmat - temp_kermat(j, k) * temp_X_diff;
          temp_part2 = temp_part2 + temp_kermat(j, k) * (X_temp.head(p) - Xbar_temp.row(j).head(p).t());
        }
      }
      //cout << "end" << endl;
    }

    //cout << temp_part1 << endl;
    //cout << temp_part2 << endl;

    res = res + temp_part2;
    //res_temp_ind  = temp_part1;
  }

  return res ;
}


// [[Rcpp::export]]
arma::mat EEZV_c(const arma::rowvec &gamma,
                const arma::rowvec &theta,
                Rcpp::ListOf<NumericMatrix> &kerMat,
                Rcpp::ListOf<NumericVector> &meas_times,
                Rcpp::ListOf<NumericMatrix> &covariates,
                Rcpp::ListOf<NumericMatrix> &Xbar,
                Rcpp::ListOf<NumericVector> &XXtbar,
                Rcpp::ListOf<NumericVector> &response,
                const arma::mat &Hmat,
                const arma::mat &Amat,
                const arma::vec &censor,
                const unsigned int &n,
                const unsigned int &p)
{

  unsigned int i = 0;
  arma::mat temp_kermat, temp_cov, Xbar_temp, temp_response;
  arma::vec X_temp, temp_part1, temp_part2, temp_part3 = vec(p + 1);
  arma::mat temp_X_diff = mat(p + 1, p + 1);
  arma::cube XXtbar_temp;

  arma::rowvec expgammaZ;

  arma::mat res_temp_ind = vec(p + 1);
  arma::mat Bmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat Sigmat = mat(p + 1, p + 1, arma::fill::zeros);
  arma::mat res = mat(p, p, arma::fill::zeros);
  arma::vec temp_meas_time_i, temp_meas_time_l, temp_countprocess_l;

  arma::vec censorind, temp_S0;

  arma::uvec temp_uvec;

  arma::uvec temp_counting;

  arma::mat HAi = Hmat * inv(Amat);

  unsigned int j, k, Jn, Kn;
  X_temp = vec(p + 1, arma::fill::zeros);
  for (i = 0; i < n; i++)
  {
    //cout << "i: " << i << endl;
    temp_meas_time_i = vec(meas_times[i].begin(), meas_times[i].size(), false);


    temp_cov = mat(covariates[i].begin(), covariates[i].nrow(),
                   covariates[i].ncol(), false);
    temp_kermat = mat(kerMat[i * n + i].begin(), kerMat[i * n + i].nrow(),
                      kerMat[i * n + i].ncol(), false);

    Xbar_temp = mat(Xbar[i].begin(), Xbar[i].nrow(),
                    Xbar[i].ncol(), false);

    Jn = temp_meas_time_i.size();
    Kn = temp_cov.n_cols;

    //XXtbar_temp = arma::cube(XXtbar[i].begin(), p + 1, p + 1, Jn, false);

    temp_part1 = vec(p, arma::fill::zeros);
    temp_part2 = vec(p, arma::fill::zeros);

    for (j = 0; j < Jn; j++){
      for (k = 0; k < Kn; k++){
        if (abs(temp_kermat(j, k)) > 1e-10){
          //Bmat = Bmat - temp_kermat(j, k) * temp_X_diff;
          temp_part1 = temp_kermat(j, k) * (temp_cov.col(k) - Xbar_temp.row(j).head(p).t());
          temp_part2 = temp_part2 + temp_part1 * temp_part1.t();
        }
      }
    }

  res = res + temp_part2;
    //res_temp_ind  = temp_part1;
  }

  return res;
}


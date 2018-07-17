#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
#include "util.hpp"


using namespace arma;
using namespace Rcpp;
using namespace std;

// [[Rcpp::export]]
NumericMatrix outermin_C(NumericVector x,NumericVector y) {
  int i,j;
  int nx =x.size();
  int ny = y.size();
  NumericMatrix res(nx, ny);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      res(i,j)=x[i]-y[j];
    }
  }
  return res;
}

// Epanechnikov Kernel function
// [[Rcpp::export]]
double epanker_C(double u,double h){
  double temp = u/h;
  double res = std::fabs(temp)<1 ? 0.75*(1.0-temp*temp)/h : 0;
  return res;
} 

// K(T_{ij}-R_{lu})
// [[Rcpp::export]]
NumericMatrix outerker_C(NumericVector x,
                         NumericVector y,
                         double& h) {
  int i,j;
  int nx =x.size();
  int ny = y.size();
  NumericMatrix res(nx, ny);
  for(i=0;i<nx;i++){
    for(j=0;j<ny;j++){
      //res(i,j)=(*kernel)(x[i]-y[j],h);
      res(i,j)=epanker_C(x[i]-y[j],h);
    }
  }
  return res;
}



// Generate Kernel matrix for K(T_{ij}-R_{lu}) for all i,j,l and u.
// [[Rcpp::export]]
Rcpp::List kerMatgen_C(Rcpp::ListOf<NumericVector>& meas_times,
                 Rcpp::ListOf<NumericVector>& obscov_times,
                 double& h) {
  int n = meas_times.size();
  int i,j;
  List res(n*n);
  for(i=0;i<n;i++){
    for(j=0;j<n;j++){
      res(i*n+j) = outerker_C(meas_times[i],obscov_times[j],h);
    }
  }
  return res;
}

/*** R
#timesTwo(42)

*/

#ifndef _NUMERIC_HPP
#define _NUMERIC_HPP

#include <RcppArmadillo.h>
//[[Rcpp::depends(RcppArmadillo)]]
//#include "optim.hpp"
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;


template <class T>
arma::vec broyden(const arma::vec &x_init, T &F, bool check);
template <class T>
void lnsrch(arma::vec &x0, double &f0, arma::vec &x1, double &f1, arma::vec &F1,
            arma::vec &deltaf, arma::vec &p, double &stpmax, bool &check, T &F);




template <class T>
arma::vec broyden(const arma::vec &x_init, T &F, bool check){

  // Numeric recipes parameters
  double EPS=std::numeric_limits<double>::epsilon();
  //Rcpp::Rcout << EPS << std::endl;
  double TOLF=1.0e-16, TOLX=EPS, STPMAX=10000.0, TOLMIN=1.0e-16;
  int MAXITS = 1000;
  check = false;
  unsigned int n = x_init.size(); // dimension, perhaps using p is better?
  arma::vec x1 = x_init;
  arma::vec F1 = F(x1);
  double f1 = arma::dot(F1,F1)/2;
  double f0 = 0;

  // check if initial guess is a root
  if(arma::max(abs(F1))<0.01*TOLF){
    check = false;
    return x1;
  }

  // Calculate stpmax for line searches
  double stpmax = STPMAX * std::max(arma::norm(x1,2),n*1.0);
  arma::vec s,x0,F0,w,p, deltaf = arma::zeros<vec>(n);
  bool restart = true;
  double max_temp,temp;
  arma::mat B =  arma::eye(n,n);
  unsigned int its=1, i=0;

  for(its=1;its<=MAXITS;its++){
    //Rcpp::Rcout << "Broyden method iteration: " << its << "=======================" << std::endl;

    if(restart){

      B = fdjac(x1,F1,F);
      //B = arma::eye(n,n);
      if(std::abs(arma::rcond(B))<1e-8) std::cout << "The Jacobian matrix is singular" << std::endl;
    } else {
      s = x1-x0;
      w = F1 -F0 - B *s;
      // Do not update with noisy componenets of w
      w.elem(find(arma::abs(w) < EPS*(arma::abs(F1)+arma::abs(F0)))).zeros();
      // Update B, the Jacobian
      if(arma::any(w)) B = B + w * s.t() /arma::dot(s,s);
      if(std::abs(arma::rcond(B))<1e-8) {
        std::cout << "The Jacobian matrix is singular" << std::endl;
        //restart = true;
      }
    }
    deltaf = B * F1;
    p = -1.0 * solve(B,F1);
    x0 = x1;
    F0 = F1;
    f0 = f1;
    //Rcpp::Rcout << "x0: " << x0.t() << "f0: " <<f0 << std::endl;
    //Rcpp::Rcout << "p: " << p.t() << std::endl;
    //B.print("B: ");
    // retuen new x 1 , f 1 , and F 1
    lnsrch(x0, f0, x1, f1, F1, deltaf, p, stpmax, check, F);
    //Rcpp::Rcout << "Done search: x1: " " << x1 << std::endl;
    if(arma::max(arma::abs(F1))<TOLF){
      //std::cout << "code 1" << std::endl;
      check = false;
      return x1;
    }
    if(check){ // True if line seach failed to find a new x1
      if(restart){
        //std::cout << "Failure; already tried reinitialzing the Jabobian" << std::endl;
        return x1;
      } else {
        max_temp = 0;
        for(i=0;i<n;i++){
          temp = std::abs(deltaf(i))/std::max(std::abs(x1(i)),1.0);
          if(temp >max_temp) max_temp = temp;
        }
        if(max_temp<TOLMIN){
          //std::cout << "code 2" << std::endl;
          return x1;
        } else {
          restart = true;
        }
      }
    } else {// Successful step; will use Broyden update for next step
      //std::cout << "Successful step" << std::endl;

      restart = false;
      // Test for convergence on \delta x
      max_temp = 0;
      for(i=0;i<n;i++){
        temp = std::abs(x1(i)-x0(i)) / std::max(std::abs(x1[i]),1.0);
        if(temp >max_temp) max_temp = temp;
      }
      if(max_temp<TOLX){
        //std::cout << "code 0" << std::endl;
        return x1;
      }
    }
  }
  std::cout << "MAXITS exceeded in broyden" << std::endl;
  return x1;
}

template <class T>
arma::mat fdjac(arma::vec &x, arma::vec &Fval, T &F){
  unsigned int n = x.size();
  arma::mat Jacobian = mat(n,n);
  unsigned int i;
  arma::vec xh = x;
  arma::vec Fh = arma::zeros<arma::vec>(n);
  double h;
  //h = 1.0e-8;
  for(i=0;i<n;i++){
    h = 1.0e-8 * abs(x(i));
    if (h==0.0) h = 1.0e-8;
    xh(i) = x(i) + h;
    h = xh(i) - x(i);
    Fh = F(xh);
    xh(i) = x(i);
    Jacobian.col(i) =  (Fh-Fval)/h;
  }
  return Jacobian;
}


template <class T>
void lnsrch(arma::vec &x0, double &f0, arma::vec &x1, double &f1, arma::vec &F1,
            arma::vec &deltaf, arma::vec &p, double &stpmax, bool &check, T &F){
  double ALF=1.0e-4, TOLX=std::numeric_limits<double>::epsilon();
  check = false;
  double temp = arma::norm(p,2);
  if(temp>stpmax){
    p=p/temp * stpmax;  // scale if attempted step is too big
    //Rcpp::Rcout << "Rescale p: " << p.t() << std::endl;
  }
  double gprime0 = arma::dot(deltaf, p);
  // Compute \lambda_{min}

  double lambda_min = TOLX/std::min(arma::max(arma::abs(p)),arma::max(arma::abs(p)/arma::abs(x0)));
  double lambda=1,lambda2 = 1;
  unsigned iter = 1;
  double rhs1=0, rhs2=0, f2=0, a=0, b=0, disc=0,lambda_temp=1;
  //Rcpp::Rcout << "Line search begin -------------------------" << std::endl;
  while(iter <  100){
    x1 = x0 + lambda * p;
    F1 = F(x1);
    f1 = arma::dot(F1,F1)/2;
    // Rcpp::Rcout << "Line search iter: " << iter <<endl;
    // Rcpp::Rcout << "x1: " << x1.t() ;
    // Rcpp::Rcout << "f1: " << f1 << ", f0: " << f0 << " ; " ;
    // Rcpp::Rcout << "lam: " << lambda << "gp: " << gprime0  <<endl;
    if(lambda < lambda_min){
      x1 = x0;
      check = true;
      return;
    } else {
      if(f1 <= f0+ALF*lambda*gprime0){
        return;
      } else {
        if(lambda==1){
          lambda_temp = - gprime0/(2*(f1-f0-gprime0));
        } else {
          rhs1 = f1 - f0 - lambda*gprime0;
          rhs2 = f2 - f0 - lambda2*gprime0;
          a = (rhs1/(lambda*lambda)-rhs2/(lambda2*lambda2))/(lambda-lambda2);
          b = (-lambda2*rhs1/(lambda * lambda) + lambda * rhs2 / (lambda2 * lambda2))/(lambda-lambda2);
          if(a==0){
            lambda_temp = -gprime0/(2*b);
          } else {
            disc = b*b-3*a*gprime0;
            if(disc < 0){
              lambda_temp = lambda/2;
            } else {
              if(b <= 0){
                lambda_temp = (-b+sqrt(disc))/(3*a);
              } else {
                lambda_temp = -gprime0/(b+sqrt(disc));
              }
            }
          }
          if(lambda_temp>0.5 * lambda) lambda_temp = 0.5*lambda;
        }
      }
    }
    lambda2 = lambda;
    f2 = f1;
    lambda = std::max(lambda_temp,0.1*lambda);
    iter++;
  }
}

#endif

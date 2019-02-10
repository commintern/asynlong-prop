#include <RcppArmadillo.h>

using namespace Rcpp;
using namespace std;
// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export]]
int check(arma::vec x, double d) {
  arma::uvec res;
  int res0;
  res = arma::find(x>=d,1);

  cout << res <<endl;
  return res.n_elem == 0 ? x.n_elem : res(0);
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
check(c(1.1,2.1,3.1),3)
*/

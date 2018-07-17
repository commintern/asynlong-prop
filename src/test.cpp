#include <RcppArmadillo.h>
#include <RcppEigen.h>

//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppArmadillo)]]
#include <iostream>
using namespace arma;
using namespace Rcpp;
using namespace std;

using Eigen::Map;                       // 'maps' rather than copies
using Eigen::MatrixXd;                  // variable size matrix, double precision
using Eigen::VectorXd;                  // variable size vector, double precision
using Eigen::SelfAdjointEigenSolver;    // one of the eigenvalue solvers


// Test functions
// [[Rcpp::export]]
arma::uvec foo(arma::vec z) {
  arma::uvec res = z > 0.5;
  return res;
}

// // [[Rcpp::export]]
// NumericMatrix outerprotest1(NumericVector x,NumericVector y) {
//   int i,j;
//   int nx =x.size();
//   int ny = y.size();
//   NumericMatrix res(nx, ny);
//   for(i=0;i<nx;i++){
//     for(j=0;j<ny;j++){
//       res(i,j)=x[i]*x[j];
//     }
//   }
//   return res;
// }
// 
// // [[Rcpp::export]]
// arma::mat outerprotest2(arma::vec x, arma::vec y) {
//   int nx =x.size();
//   int ny = y.size();
//   arma::mat res = mat(nx, ny);
//   res = x*y.t();
//   return res;
// }
// 
// 
// 
// // [[Rcpp::export]]
// MatrixXd outerprotest3(Map<VectorXd> X,Map<VectorXd> Y) {
//   return X * Y.adjoint();
// }
// 
// // [[Rcpp::export]]
// arma::mat mattest(Rcpp::ListOf<NumericMatrix>& kerMat) {
//   //NumericMatrix temp(kerMat[1]);
//   //double* temp = );
//   mat temp_kermat;
//   for(int i=0;i<10;i++){
//     temp_kermat =mat(kerMat[i].begin(), kerMat[i].nrow(), kerMat[i].ncol(), false);
//   }
//   return temp_kermat;
// }
// 
// // [[Rcpp::export]]
// arma::mat mattest1(const Rcpp::ListOf<NumericMatrix>& kerMat) {
//   mat temp_kermat;
//   for(int i=0;i<10;i++){
//     temp_kermat = Rcpp::as<arma::mat>(kerMat[i]);
//   }
//   
//   return temp_kermat;
// }
// 
// // [[Rcpp::export]]
// arma::uvec compa(const arma::vec& x, const double& va){
//   return x >= va;
// }
// 
// 
// // [[Rcpp::export]]
// NumericVector tran1(NumericVector x){
//   int n = x.size();
//   int i;
//   for(i=0;i<n;i++){
//     x[i] = pow(x[i],2);
//   }
//   return x;
// }
// 
// // [[Rcpp::export]]
// arma::vec tran2(arma::vec x){
//   x.transform([](double val){ return pow(val,2); });
//   return x;
// }


// // [[Rcpp::export]]
// arma::field<arma::mat> fieldgen(Rcpp::ListOf<NumericVector>& A,
//                         Rcpp::ListOf<NumericVector>& B) { 
//   std::size_t An = A.size();
//   std::size_t Bn = B.size();
//   unsigned int i,j;
//   arma::field<arma::mat> F(An,Bn);
//   for(i=0;i<An;i++){
//     for(j=0;j<Bn;j++){
//       F(i,j) = Rcpp::as<arma::vec>(A[i]) * Rcpp::as<arma::vec>(B[j]).t();
//     }
//   }
//   return F;
// }

// // [[Rcpp::export]]
// arma::field<arma::vec> fieldout(arma::field<arma::vec>& A) { 
//   std::size_t An = A.size();
//   
//   unsigned int i;
//   arma::field<arma::vec> F(An);
//   for(i=0;i<An;i++){
//     F(i) = A(i) *2;
//   }
//   return F;
// }

// // [[Rcpp::export]]
// arma::field<arma::field<arma::vec>> fielde(arma::field<arma::field<arma::vec>>& A) { 
//   cout << A(0,0) << endl;
//   return A;
// }



// // [[Rcpp::export]]
// arma::field<arma::vec> temptest1( Rcpp::ListOf<NumericVector>& meas_times) {
//   unsigned int i;
//   arma::vec temp_meas_time;
//   int n = meas_times.size();
//   arma::field<arma::vec> res(n);
//   for (i = 0; i < n; i++) {
//     temp_meas_time = vec(meas_times[i].begin(),meas_times[i].size(),false);
//     res(i)  = temp_meas_time*2;
//   }
//   return res;
// }

// // [[Rcpp::export]]
// arma::uvec countprofun_C(arma::vec counttime, arma::vec externalTime){
//   arma::uword i,j;
//   arma::uword Nc = counttime.size();
//   arma::uword Nt = externalTime.size();
//   arma::uvec res(Nt);
//   j = 0;
//   i = 0;
//   while((i<Nt)&&(j<Nc)){
//     if(j>=Nc-1){
//       res(i) = Nc;
//       i++;
//     } else {
//       if(counttime(j)<externalTime(i)){
//         if(counttime(j+1)>=externalTime(i)){
//           res(i) = j+1;
//           ++i;
//         } else {
//           ++j;
//         }
//       } else {
//         res(i) = j;
//         ++i;
//       }
//     }
//     
//   }
//   return res;
// }


// // Compare field and cube
// // [[Rcpp::export]]
// arma::field<arma::vec> filrd(arma::vec vv, int nrow, int ncol){
//   int i,j;
//   int p = vv.size();
//   arma::field<arma::vec> res(nrow,ncol);
//   arma::vec temp = vec(p+1);
//   for(i=0;i<nrow;i++){
//     for(j=0;j<ncol;j++){
//       temp.subvec(0,p-1)=vv;
//       temp(p) = 1;
//       res(i,j)=temp;
//     }
//   }
//   
//   for(i=0;i<nrow;i++){
//     for(j=0;j<ncol;j++){
//       temp=res(i,j);
//     }
//   }
//   return res;
// }
// 
// // [[Rcpp::export]]
// arma::cube cuberd(arma::vec vv, int nrow, int ncol){
//   int i,j;
//   int p = vv.size();
//   arma::cube res(nrow,ncol,p+1);
//   arma::vec temp = vec(p+1);
//   for(i=0;i<nrow;i++){
//     for(j=0;j<ncol;j++){
//       res.subcube( i, j, 0, i, j, p-1 )=vv;
//       res.subcube( i, j, p, i, j, p )=1;
//     }
//   }
//   
//   for(i=0;i<nrow;i++){
//     for(j=0;j<ncol;j++){
//       temp=res.tube(i,j);
//     }
//   }
//   return res;
// }

// [[Rcpp::export]]
arma::cube Xgen_C(const arma::mat& covMat, const arma::vec& countprocess, int p){
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

/*** R
#testres <- microbenchmark::microbenchmark(outerprotest1(1:100,1:100),outerprotest2(1:100,1:100),outerprotest3(as.numeric(1:100),as.numeric(1:100)))
#plot(testres)

#tm <- lapply(1:10, function(x) matrix(runif(30*30),30))
#tm
#mattest(tm)
#mattest1(tm)

#tr <- microbenchmark::microbenchmark(mattest(tm),mattest1(tm),times=1000L)  
#plot(tr)
#tt <- runif(10)
#compa(tt,tt[1])
#microbenchmark::microbenchmark(tran1(1:200),tran2(1:200),times=100L)  
#t1 <- replicate(2,runif(2),simplify = F)
#t2 <- replicate(2,runif(2),simplify = F)
#test <- fieldgen(t1,t2)
#t1
#fieldout(t1)
#microbenchmark::microbenchmark(fieldout(meas_obs_list),temptest1(meas_obs_list))
#microbenchmark::microbenchmark(cuberd(1:4,6,7),filrd(1:5,6,7))
*/

#ifndef LONGEST_H
#define LONGEST_H

#include <RcppArmadillo.h>
arma::vec countprofun_C(const arma::vec& counttime, const arma::vec& externalTime);

arma::cube Xgen_C(const arma::mat& covMat, const arma::vec& countprocess, const unsigned int& p);
#endif

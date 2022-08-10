#ifndef cellMCD_H
#define cellMCD_H

#ifndef ARMA_DONT_PRINT_ERRORS
#define ARMA_DONT_PRINT_ERRORS
#endif

#ifndef  ARMA_USE_CXX11
#define ARMA_USE_CXX11
#endif

#include "RcppArmadillo.h" 

arma::mat subinverse_cpp(const arma::mat& Sigma,
                         const arma::mat& Sigmai,
                         arma::uvec indx);

void uniqueRows(arma::umat W,
                std::unordered_map<std::string, arma::uvec>& Wmap);

arma::vec Deltacalc_cpp(const arma::mat & X,
                        const arma::umat & W,
                        const arma::mat & Sigma,
                        const arma::mat & Sigmai,
                        const arma::vec & mu,
                        arma::uword j);


#endif

#include <armadillo>
#include <functional>

#ifndef MCMC_H
#define MCMC_H

 
void run_mcmc( std::function<double(arma::vec&)> f, arma::vec& X, long int N,
	       std::function<void(arma::vec&)> obs);
#endif /* MCMC_H */

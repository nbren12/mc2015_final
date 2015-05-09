#include <armadillo>

#ifndef MCMC_H
#define MCMC_H

template<typename func, typename observer>
  void run_mcmc(func& f, arma::vec& X, long int N, observer obs);

#endif /* MCMC_H */

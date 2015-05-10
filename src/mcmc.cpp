#include <armadillo>
#include <iostream>
#include <gsl/gsl_rng.h>
#include <gsl/gsl_randist.h>

#include "mcmc.h"

#define NUMP 3

gsl_rng * r;

double stdev[]= {.5, .1, 1.0};

using namespace arma;
using namespace std;

void proposal_rand(vec& X, vec& Y){
  int i;
  for (i = 0; i < NUMP; i++) {
    Y(i) = X(i) + gsl_ran_gaussian(r, stdev[i]);
  }
}

// Proposal step generate gaussian
double proposal_pdf(vec& X, vec& Y){
  int i;
  double p = 1.0;
  for (i = 0; i < NUMP; i++) {
    p *= gsl_ran_gaussian_pdf(Y(i)-X(i), stdev[i]);
  }

  return p;
}

// Metropolis accept function
template<typename func> int A(func& f, vec& X, vec& Y ){
  double accept;
  accept = f(Y) * proposal_pdf(Y, X) /
    (f(X) * proposal_pdf(X, Y));

  if (accept >= 1){
    return 1;
  } else {
    return gsl_ran_bernoulli(r, accept);
  }
}

// MCMC runner
void run_mcmc( std::function<double(vec&)> f, vec& X, long int N,
	       std::function<void(vec&)> obs){
  int i;
  r = gsl_rng_alloc(gsl_rng_default);
  // Copy X into new dvec Y
  vec Y = X;

  for (i = 0; i < N; i++) {
    cout << "Sample "<< i << " of " << N << endl;
    proposal_rand(X, Y);
    if (A(f, X, Y))
      X = Y;

    // Observe X
    obs(X);
  }

  gsl_rng_free(r);
}

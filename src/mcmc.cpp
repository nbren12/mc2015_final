#include <armadillo>
#define gaussian_length_scale 1.0

void proposal_rand(vec& X, vec& Y){
  int d = X.n_rows;
  int i;
  for (i = 0; i < d; i++) {
    Y(i) = X(i) + gsl_ran_gaussian(r, gaussian_length_scale);
  }
}

// Proposal step generate gaussian
double proposal_pdf(vec& X, vec& Y){
  int d = X.n_rows;
  int i;

  double p = 1.0;
  for (i = 0; i < d; i++) {
    p *= gsl_ran_gaussian_pdf(Y(i)-X(i), gaussian_length_scale);
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
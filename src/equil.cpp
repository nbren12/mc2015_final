#include <armadillo>
#include <functional>
#include <cmath>
#include "integrate.h"


double equil(vec& theta, mat& U, vec& t, const double beta){
  double prob = 0.0;
  
  // Work array
  vec y(U.n_rows);
  int i;
  for (i = 0; i < t.n_elem-1; i++) {
    y = evolve(t(i), t(i+1), U.col(i), theta);
    prob += norm(y-U.col(i+1), 2);
  }

  prob = exp(-beta * prob/ (t.n_elem -1));

  return prob;
  
}



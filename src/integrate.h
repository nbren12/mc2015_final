#include <functional>
#include <armadillo>


#ifndef INTEGRATE_H
#define INTEGRATE_H


using namespace arma;

arma::vec evolve(double tStart, double tEnd, arma::vec y, arma::vec& theta,
		 double dtMax = .001);



mat run_model(vec& tout, vec & y, double dtMax = .001);

#endif /* INTEGRATE_H */



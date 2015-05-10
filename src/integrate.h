#include <functional>
#include <armadillo>


#ifndef INTEGRATE_H
#define INTEGRATE_H


using namespace arma;

arma::vec evolve(double tStart, double tEnd, arma::vec y, arma::vec& theta,
		 double dtMax = .01);



mat run_model(vec& tout, vec & y, double dtMax=.001,
	      double sigma=10.0, double beta=8.0/3.0, double r=28.0, double epsilon=0.0);

#endif /* INTEGRATE_H */



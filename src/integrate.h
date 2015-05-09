#include <functional>
#include <armadillo>

using namespace arma;

arma::vec evolve(double tStart, double tEnd, arma::vec y, arma::vec& theta,
		 double dtMax);



mat run_model(vec& tout, vec & y, double dtMax = .001);

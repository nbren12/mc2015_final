#include <functional>
#include <armadillo>

using namespace arma;

arma::vec evolve(double tStart, double tEnd, arma::vec y, arma::vec& theta,
		 double dtMax);


std::function<double(vec&)> mk_equil_functor(mat& U, vec& t);

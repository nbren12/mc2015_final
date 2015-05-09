#include <armadillo>
#include "mcmc.h"
#include "equil.h"
#include "integrate.h"

using namespace arma;

int main(int argc, char *argv[])
{
  double tau  = 1.0;   // Sampling interval
  double tEnd = 100.0; // Sampling length

  vec tout = linspace(0.0, tEnd, (int) tEnd /tau);	// Output times
  vec y("-7.4185 -12.4638 15.6519");			// initial condition

  cout << "Running Truth Model" << endl;
  auto output = run_model(tout, y);			// get output data

  vec theta; // Parameters
  theta << 10.0 << 8.0/3.0 << 28.0 << 0.0;

  cout << equil(theta, output, tout, 1.0);

  return 0;
}

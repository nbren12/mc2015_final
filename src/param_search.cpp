#include <functional>
#include <armadillo>
#include "mcmc.h"
#include "equil.h"
#include "integrate.h"

using namespace arma;
using namespace std;
using namespace std::placeholders;

int test_equil_distr()
{
  double tau  = 1.0;   // Sampling interval
  double tEnd = 100.0; // Sampling length

  vec tout = linspace(0.0, tEnd, (int) tEnd /tau);	// Output times
  vec y("-7.4185 -12.4638 15.6519");			// initial condition

  cout << "Running Truth Model" << endl;
  auto output = run_model(tout, y);			// get output data

  vec theta;					// Parameters
  theta << 10.0 << 8.0/3.0 << 28.0 << 0.0;	// True values
  theta(0) = 9.5;				// Make sigma different
  cout << equil(theta, output, tout, .01);      // Print equil output

  return 0;
}

int main(int argc, char *argv[])
{
  long int N  = 10000;   // Number of samples
  double noise = 1.0;  // Imperfect model
  double tau  = 10.0;   // Sampling interval
  double beta = .04;   // Inverse temperature
  double tEnd = 100.0; // Sampling length

  vec tout = linspace(0.0, tEnd, (int) tEnd /tau);	// Output times
  vec y("-7.4185 -12.4638 15.6519");			// initial condition

  cout << "Running Truth Model" << endl;
  auto output = run_model(tout, y);			// get output data

  vec theta;					// Parameters
  theta << 10.0 << 8.0/3.0 << 28.0 << noise;	// True values
  theta(0) = 10.5;				// Make sigma different

  auto f = bind(equil, _1, output, tout, beta);  // equildist for theta f(theta) 

  // Setup output file
  ofstream outfile;
  outfile.open("samples.txt");
    
  run_mcmc(f, theta, N, [&](vec& x){
      int i;
      for (i = 0; i < x.n_elem; i++) {
	outfile << x[i];
	if (i < x.n_elem -1) outfile << " ";
      }
      outfile << endl;
    });
  outfile.close();
  return 0;
}

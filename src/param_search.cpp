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

void estimate_parameters(long int N, double noise, double beta, double tau, int nt) {

  // Make a vector of output times at regular intervals set by tau
  vec tout(nt+1);	
  int i;
  for (i = 0; i < nt +1; i++) {
    tout(i) = i *tau;
  }

  vec y("-7.4185 -12.4638 15.6519");		// initial condition of Lorenz63 system
  auto output = run_model(tout, y);		// Get truth model samples for times in tout

  vec theta;					// Vector of parameters
  theta << 10.0 << 8.0/3.0 << 28.0 << noise;	// True values of parameters
  theta << 7.0 << 5.0/3.0 << 20.0 << noise;	// init cond of parameters

  auto f = bind(equil, _1, output, tout, beta); // equildist for theta f(theta) 

  // Setup output file
  ofstream outfile;
  outfile.open("samples.txt");
  outfile << "sigma,beta,r,noise" << endl;

  // Run the sampler
  run_mcmc(f, theta, N, [&](vec& x){
      int i;
      for (i = 0; i < x.n_elem; i++) {
	outfile << x[i];
	if (i < x.n_elem -1) outfile << ",";
      }
      outfile << endl;
    });
  outfile.close();
}

int main(int argc, char *argv[])
{
  int offset = 1;
  double tau = atof(argv[++offset]);
  int nt = atol(argv[++offset]);

  if  (strcmp(argv[1], "sample") == 0){
    double noise = atof(argv[++offset]);
    double beta = atof(argv[++offset]);  // Gibbs distribution parameter
    long int N = atol(argv[++offset]);
    estimate_parameters(N,  noise,  beta,  tau, nt);

  } else if  (strcmp(argv[1], "run") == 0){

    double sigma = atof(argv[++offset]);
    double beta = atof(argv[++offset]);  // Parameter in Lorenz system
    double r = atof(argv[++offset]);
    double epsilon = atof(argv[++offset]);

    vec tout = linspace(0, nt * tau, nt+1);             // Output times
    vec y("-7.4185 -12.4638 15.6519");			// initial condition
    auto output = run_model(tout , y, .001, sigma, beta, r, epsilon);
    output.save("output.txt", raw_ascii);
  }
  return 0;
}

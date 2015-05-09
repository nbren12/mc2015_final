#include <iostream>
#include <stdlib.h>
#include <armadillo>
#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

using namespace std;
using namespace arma;


/* extern SdeCoeff sdeA, sdeB; */
/* typedef  int (* SdeCoeff)(double, const double[], double[], void * ); */

struct WorkStruct {
  const int n;
  vec a, b, c;
  double  sigma, beta, rho, epsilon;
  gsl_rng *rng;

  // Constructor
  WorkStruct(int n, double sigma=10.0, double beta=8.0/3.0, double r=28.0,
	     double epsilon=0.0) :
    n(n), a(vec(n)), b(vec(n)), c(vec(n)),
    sigma(10.0), beta(beta), rho(r), epsilon(epsilon) {
    rng = gsl_rng_alloc(gsl_rng_default);
  }

  // Destructor
  ~WorkStruct(){
    free(rng);
  }

  // time stepping routines
  void eulerStep(double t, double dt, vec& y);
  void evolve(double tStart, double tEnd, vec& y,
	      double dtMax);

  // Definition of SDE (can be refactored out)
  int sdeA(double t,  double y[], double dydt[]) {

    dydt[0] = sigma * (y[1] - y[0]);
    dydt[1] = y[0] * ( rho - y[2]) - y[1];
    dydt[2] = y[0] * y[1] - beta * y[2];

    return 0;
  }

  /* Noise term for  the Lorenz63 system */
  int sdeB(double t, const double y[], double b[]) {

    int i;
    for (i = 0; i < n; i++) {
      b[i] = epsilon;
    }

    return 0;
  }
  
};

void WorkStruct::eulerStep(double t, double dt, vec& y){
  sdeA(t, y.memptr(), a.memptr());

  // Deterministic predictor corrector
  b = y + dt * a;
  sdeA(t, b.memptr(), c.memptr());
  y = y  + dt/2.0 * (a + c);


  // Stochastic step forward euler
  sdeB(t, y.memptr(), a.memptr());
  int i;
  double ys;
  for (i = 0; i < n; i++) {
    y[i] += gsl_ran_gaussian(rng, sqrt(dt)*b(i));
  }
}

void WorkStruct::evolve(double tStart, double tEnd, vec& y, double dtMax){

  double dt, t;
  t = tStart;
  dt = dtMax;
  int killLoop = 0;
  
  /* Take timesteps of size dtMax until tEnd is reached */
  while (1){

    if (dt > 0) {
      eulerStep(t, dt, y);
      t += dt;
      cout << t << endl;
    }
    if (killLoop) break;

    if ( t > tEnd - dt) {
      dt = tEnd - t;
      killLoop = 1;
    }
  }
}

template<typename stream, typename T> void save(stream& fout, T & y){
  int i;
  for (i = 0; i < y.n_elem; i++) {
    fout << y(i);
    if ( i < y.n_elem - 1) fout << " ";
  }
  fout << endl;
}


int test_eulerStep()
{
   // initialize working arrays
  const int n =3;

  WorkStruct work(n);

  vec y;
  y << 1 << 1 << 1;
  double t = 0;
  double dt = .1;

  // step function
  work.eulerStep(t, dt, y);

  return 0;
}



int test_evolve()
{
   // initialize working arrays
  const int n =3;
  WorkStruct work(n);
  vec y;
  y << 1 << 1 << 1;
  double tStart = 0;
  double tEnd   = 1.0;
  double dtMax = .01;

  // Evolve function from start to finish
  work.evolve(tStart, tEnd, y, dtMax);

  return 0;
}

int test_output()
{
   // initialize working arrays
  const int n =3;
  WorkStruct work(n);
  vec y;
  y << 1 << 1 << 1;
  double tStart = 0;
  double tEnd   = 1.0;
  double dtMax = .01;

  // Output times
  vec tout = linspace<vec>(tStart, tEnd, 100);

  // Output file
  ofstream fout;
  fout.open("output.txt");
  save(fout, y);

  int i;
  for (i = 1; i < tout.n_elem; i++) {
    work.evolve(tout(i-1), tout(i), y, dtMax);
    save(fout, y);
  }
  fout.close();
  return 0;
}


int main(int argc, char *argv[])
{
  // test_eulerStep();
  // test_evolve();
  test_output();
  return 0;
}

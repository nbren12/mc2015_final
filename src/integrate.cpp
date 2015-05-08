#include <iostream>
#include <stdlib.h>
#include <armadillo>
#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "header.h"

using namespace std;
using namespace arma;


/* extern SdeCoeff sdeA, sdeB; */
/* typedef  int (* SdeCoeff)(double, const double[], double[], void * ); */

struct WorkStruct {
  const int n;
  vec a, b, c;

  gsl_rng *rng;

  WorkStruct(int n) : n(n), a(vec(n)), b(vec(n)), c(vec(n)) {
    rng = gsl_rng_alloc(gsl_rng_default);
  }
  ~WorkStruct(){
    free(rng);
  }

  void eulerStep(double t, double dt, vec y);
  void evolve(double tStart, double tEnd, vec& y,
	      double dtMax);
};

void WorkStruct::eulerStep(double t, double dt, vec y){
  sdeA(t, y.memptr(), a.memptr(), NULL);

  // Deterministic predictor corrector
  b = y + dt * a;
  sdeA(t, b.memptr(), c.memptr(), NULL);
  y = y  + dt/2.0 * (a + c);


  // Stochastic step forward euler
  sdeB(t, y.memptr(), a.memptr(), NULL);
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


int main(int argc, char *argv[])
{
  test_eulerStep();
  test_evolve();
  return 0;
}

#include <iostream>
#include <stdlib.h>
#include <math.h>
#include "gsl/gsl_rng.h"
#include "gsl/gsl_randist.h"

#include "header.h"

using namespace std;


/* extern SdeCoeff sdeA, sdeB; */
/* typedef  int (* SdeCoeff)(double, const double[], double[], void * ); */

struct WorkStruct {
  int n;
  double *a, *b;

  gsl_rng *rng;

};

struct WorkStruct work;

void initializeWorkStruct(int n){
  work.n = n;
  work.a = (double*) malloc(n * sizeof(double));
  work.b = (double*) malloc(n * sizeof(double));
  work.rng = gsl_rng_alloc(gsl_rng_default);
}

void destroyWorkStruct(){
  free(work.a);
  free(work.b);
  free(work.rng);
}

void eulerStep(double t, double dt, double y[]){
  sdeA(t, y, work.a, NULL);
  sdeB(t, y, work.b, NULL);

  int i;
  for (i = 0; i < work.n; i++) {
    y[i] += dt * work.a[i] +
            gsl_ran_gaussian(work.rng, sqrt(dt)*work.b[i]);
  }
}

void evolve(double tStart, double tEnd, double y[], double dtMax){

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
  initializeWorkStruct(n);

  double y[n] = { 1, 1, 1};
  double t = 0;
  double dt = .1;

  // step function
  eulerStep(t, dt, y);

  return 0;
}



int test_evolve()
{
   // initialize working arrays
  const int n =3;
  initializeWorkStruct(n);

  double y[n] = { 1, 1, 1};
  double tStart = 0;
  double tEnd   = 1.0;
  double dtMax = .01;

  // Evolve function from start to finish
  evolve(tStart, tEnd, y, dtMax);

  return 0;
}


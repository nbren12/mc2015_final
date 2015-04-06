/* Define constants and parameters */
#define N 3
#define sigma 10.0
#define beta 8.0/3.0
#define rho 28.0
#define epsilon 0.5



/* The deterministic Lorenz63 system from wikipedia */
int sdeA(double t, const double y[], double dydt[], void * params) {

  dydt[0] = sigma * (y[1] - y[0]);
  dydt[1] = y[0] * ( rho - y[2]) - y[1];
  dydt[2] = y[0] * y[1] - beta * y[2];

  return 0;
}

/* Noise term for  the Lorenz63 system */
int sdeB(double t, const double y[], double b[], void * params) {

  int i;
  for (i = 0; i < N; i++) {
    b[i] = epsilon;
  }

  return 0;
}

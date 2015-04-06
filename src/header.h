#ifdef __cplusplus
extern "C" {
#endif // __cplusplus

  int sdeA(double t, const double y[], double dydt[], void * params);
  int sdeB(double t, const double y[], double dydt[], void * params);
  void evolve(double tStart, double tEnd, double y[], double dtMax);

  // Used to initialize and destroy working variables
  void initializeWorkStruct(int n);
  void destroyWorkStruct();

#ifdef __cplusplus
}
#endif // __cplusplus

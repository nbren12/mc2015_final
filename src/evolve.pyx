"""
A small wrapper module for the C function


import using pyximport(include_args


import pyximport
import numpy as np

from os import environ
environ['LDFLAGS'] = '-L./build -llorenz63'

pyximport.install(setup_args={'include_dirs':[np.get_include()]})
import evolve.evolvepy as evolve
"""
cimport numpy as np
import numpy as np

cdef extern:
  void evolve(double tStart, double tEnd, double y[], double dtMax);
  void initializeWorkStruct(int n);
  void destroyWorkStruct();

def evolvepy(double tStart, double tEnd, double[:] y, double dtMax=.001):
    evolve(tStart, tEnd, &y[0], dtMax)

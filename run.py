import subprocess
import numpy as np
import pandas as pd
EXE = 'build/lorenzmcmc'


def run(tau=.01, nt=1000, sigma=10.0, beta=8/3, r=28, noise=0.0):
    args = [EXE, 'run']
    
    floatargs = [tau, nt, sigma, beta, r, noise]
    args += list(map(str, floatargs))
    subprocess.call(args)
    t = np.arange(nt + 1) * tau
    return t, np.genfromtxt("output.txt")


def sample(N=1000, tau=1.0, nt=10, noise=.1, beta=.01):
    args = [EXE, 'sample']
    floatargs = [tau, nt, noise, beta, N]
    args += list(map(str, floatargs))
    subprocess.call(args, stdout=subprocess.PIPE)
    return pd.read_csv("samples.txt")

if __name__ == '__main__':
    print(sample())

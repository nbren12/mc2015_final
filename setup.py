import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize
import numpy as np

print("Compiling C library using cmake")
cwd =  os.getcwd()
src = os.path.join(cwd, 'src')

os.chdir(src)
os.system('cmake . && make')
os.chdir(cwd)

ext = Extension(src + '/evolve', [src + '/evolve.pyx'],
                include_dirs=[np.get_include()],
                library_dirs=[src,],
                libraries=['lorenz63'])

setup(
    name="My hello app",
    packages=['src',],
    ext_modules=cythonize(ext),
)

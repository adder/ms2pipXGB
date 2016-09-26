from distutils.core import setup
from Cython.Build import cythonize
import numpy

setup(
    ext_modules = cythonize("ms2pipfeatures_cython.pyx"),
    include_dirs=[numpy.get_include()]
)

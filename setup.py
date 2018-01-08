from distutils.core import setup, Extension
from Cython.Build import cythonize
import numpy 

setup(
    name="bilayer-clusters",
    version="0.1",
    packages=["bilayer_clusters",],
    description= "GROMACS Trajectory Analysis Repository",
    author="Tristan Ang",
    license='MIT License',
    ext_modules = cythonize("bilayer_clusters/*.pyx"),
    include_dirs=[numpy.get_include()]

)

#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.

#import distribute_setup
#distribute_setup.use_setuptools()
from setuptools import find_packages, setup
#from distutils.core import setup
from distutils.extension import Extension
from Cython.Distutils import build_ext
from Cython.Build import cythonize
#from setuptools import find_packages

import sys
import numpy

includes = [numpy.get_include()]
libpath = []

Win64 = sys.platform.startswith("win")

openmp = '/openmp' if Win64 else '-fopenmp'

ext_modules = cythonize("raypier/core/*.pyx",
                        language="c++" if Win64 else None,
                        include_path=[numpy.get_include()],
                        language_level = "3")

for module in ext_modules:
    module.extra_compile_args.append(openmp)
    module.extra_link_args.append(openmp)
    
with open("README.rst", 'r', encoding="utf-8") as fh:
    long_description = fh.read()

setup(
    name="raypier",
    version="0.2.3",
    packages=find_packages(),
    scripts = [], #no stand-alone application yet
    cmdclass = {'build_ext': build_ext},
    #ext_package = "raypier", ###Not needed on linux. Is it necessary on Win64?
    ext_modules = ext_modules,
    include_dirs = [numpy.get_include()],
    zip_safe = False, #Zipped eggs cause the ReadTheDocs build to fail.

    install_requires = [
        #"VTK",
        #"wxPython", #maybe can avoid an explicit dependancy on wx
        "numpy >= 1.19", 
        "traits >= 6.0",
        "mayavi >= 4.7", #TVTK is distributed as part of Mayavi
        "traitsui >= 7.0",
        "pyyaml"
        ],

    package_data = {
        "": ["*.pyx", "*.pxd"],
        "raypier": ["material_data/glass_dispersion_database.db"]
        },

    author = "Bryan Cole",
    author_email = "bryancole.cam@gmail.com",
    description = """A optical ray-tracing package, for design, optimisation and \
visualisation of mirror/lens systems.""",
    long_description = long_description,
    license = "GPL3",
    keywords = "science engineering optics ray-tracing physics",
    classifiers=["Development Status :: 4 - Beta",
                 "License :: OSI Approved :: GNU General Public License v3 (GPLv3)",
                 "Natural Language :: English",
                 "Programming Language :: Cython",
                 "Programming Language :: Python :: 3.7",
                 "Programming Language :: Python :: Implementation :: CPython",
                 "Topic :: Scientific/Engineering :: Physics",
                 "Topic :: Scientific/Engineering :: Visualization"
                 ],
    url = "https://groups.google.com/u/1/g/python-raytrace",
    project_urls = {
        "Homepage" : "https://groups.google.com/u/1/g/python-raytrace",
        "Documentation" : "https://raypier-optics.readthedocs.io/en/latest/index.html",
        "Source" : "https://github.com/bryancole/raypier_optics.git",
        "Issues" : "https://github.com/bryancole/raypier_optics/issues"
        },
    python_requires=">=3.7"
    )

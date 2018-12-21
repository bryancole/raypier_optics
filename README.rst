=========================================================
Raytrace - A Non-sequential optical ray-tracing framework
=========================================================

Raytrace is a non-sequential ray-tracing framework for modelling optical systems. Its main features are:

  # It's pretty fast. The tracing algorithms are written in Cython (i.e. C) and use efficient data-structures for handling many thousands of rays.
  # Correctly handles polarisation
  # Support for dispersion including a comprehensive library of glass-types (taken from refractive-index.info)
  # An interactive user-interface with 3D model/ray rendering through VTK. Tracing is performed "live".
  # A modest selection of optic-types including singlet and achromatic doublet lenses (with AR coating), plane-mirrors, 
    ellipsoidal mirrors, prisms
  # Experimental support for dispersion-analysis (i.e. ultra-fast optics)
  # (Currently broken) support for STEP export of model and rays, via PythonOCC library.
  
There are still a few "big" features missing:

  # A point-spread function. I'll implement this eventually...
  # Need to re-work STEP export
  
Requirements
============

Raytrace requires:

  # python-3.6 (older versions work on 2.7) 
  # numpy
  # traits / traitsui
  # Cython
  # Mayavi (for TVTK) / VTK
  # (Optionally) PythonOCC - for STEP export
  # (Optionally) Chaco / Enable - for integrating line-plots into the UI
  
Installation
============

The best way to install and run Raytrace is using a Conda environment. 

  # Install miniconda (or the full Anaconda distribution)
  # Create a fresh conda-environment using the environment.yml file included 
    in this raytrace repo.
Installation
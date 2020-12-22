Introduction to Raytrace
========================

Raytrace is a non-sequential optical ray-tracing program. It is intended as a 
design tools for modelling optical systems (cameras, imaging systems, telescopes etc.).

The main features of ray-trace are:
 - Non-sequential tracing (no need to specify the order of optical components)
 - Nice visualisation of the traced result
 - Live update to the traced result as the user adjusts the model
 - Reasonable performance (tracing algorithms runs at C-speed using Cython)
 - STEP export of models for integration with CAD design (using PythonOCC)
 - Saving / Loading models in YAML format.
 - Trace rays with full polarisation and phase information
 - Dielectric Materials with simple-coating supported, including dispersion
 - A basic library of materials (from RefractiveIndex.info)
 - 
 - Various analysis algorithms including E-field evaluation by sum-of-Gaussian-Modes, and
   dispersion calculations for ultra-fast optics applications.

At present, the primary means of using raytrace is to create your model in the
form of a python script. However, it is possible to launch an empty model and then 
add in components from the GUI / Menu.

A minimal empty model looks like::

  from raytrace.tracer import RayTraceModel
  model = RayTraceModel()
  model.configure_traits()

This opens a GUI window from which you can add model objects using the Insert menu.


Build and Installation
======================

Installing the dependencies is probably the biggest hurdle to using Raytrace. The conda
package/environment manager is far and away the easiest means of getting all the requirements
installed.

Once you have the requirements, building raytrace requires a compiler. The usual process of::

    python setup.py build
    sudo python setup.py install
    
should work [no sudo if you're on Windows]. I've never tried building on a Mac.
    
If I could figure out how the heck conda-forge worked, I'd probably use it.

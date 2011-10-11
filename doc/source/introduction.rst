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
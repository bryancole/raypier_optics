===========
Distortions
===========

The :py:mod:`raypier.distortions` module contains objects representing distortions of a given face. The Distortion objects
are intended to be used with the :py:class:`raypier.faces.DistortionFace` class (part of the General Optic framework). 
Fundamentally, any 2D function can be implemented as a Distortion. At present, on a single type is implemented. I intend 
to implement a general Zernike polynomial distortion class. Other distortion functions are easy to add.

An example of their usage:

.. literalinclude:: /../../examples/distortion_face_example.py

This example shows a high-amplitude distortion, for illustrative purposes.

.. image:: images/distortion_example.png

During the ray-tracing operation, the intersections with distorted faces are found using an iterative algorithm similar to Newton-
Ralphson. Typically, the intersection is found with 2 to 3 calls to the intercept-method of the underlying face. Distortions
are expected to be small deviations from the underlying face (maybe no more than a few wavelengths at most). If you make the amplitude of the distortion large, the under of 
iterations to converge will increase and the ray-tracing hit take a performance hit. For very large distortions, the intercept 
my fail altogether. 

One could, in principle, wrap multiple DistortionFaces over other DistortionFaces. However, I would expect the performance penalty to be 
quite severe. In this case, A better plan would be to implement a specialised DistortionList object which can sum the distortion-values
from a list of input Distortions. On my todo list ...

In python scripting, one can simply evaluate any Distortion object given some x- and y-coordinates as numpy arrays. This is 
useful for testing. For example::

    from raypier.distortions import SimpleTestZernikeJ7
    import numpy
    
    dist = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=0,1)
    
    x=y=nmupy.linspace(-10,10,500)
    X,Y = numpy.meshgrid(x,y)
    
    Z = dist.z_offset(X.ravel(),Y.ravel())
    Z.shape = X.shape #restore the 2D shape of the Z array
    
Distortions have an additional method, :py:meth:`Distortion.z_offset_and_gradient`. This returns a array of shape (N,3) 
where the input X and Y arrays have length N. The first two columns of this array contain the gradient of the 
distortion, dZ/dX and dZ/dY respectively. The third column simply contains Z. Returning both Z and it's gradient 
turns out to be useful at the C-level during tracing. I.e.::

    grad = dist.z_offset_and_gradient(X.ravel(), Y.ravel()).reshape(X.shape[0], X.shape[1],3)
    dZdX = grad[...,0]
    dZdY = grad[...,1]
    Z = grad[...,2]

Zernike Polymonial Distortions
------------------------------

More general distortions can be applied using the :py:class:`raypier.distortions.ZernikeSeries` class.

As previously, instances of this object are passed to a :py:class:`raypier.faces.DistortionFace` , along
with the base-surface to which the distortion is to be applied.
    
An example of the this class in action can be seen here:

.. literalinclude:: /../../examples/zernike_distortion_example.py

Here's what the model looks like in the UI.

.. image:: images/zernike_distortions_example.png

This example also demonstrates the use of a Constraints object to provide some UI controls for
easier adjustment of the relevant model parameters.
    

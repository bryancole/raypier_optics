============================
Gaussian Beamlet Propagation
============================

Gaussian beams are a solution to the paraxial wave equation. It turns out that the propagation of such modes
can be represented by a chief ray and a number of skew marginal rays which define the 1/e^2 edge (and divergence) of the beam.
It has been shown that such the propagation of such modes obeys the rules of geometric ray tracing (i.e. Snell's law).

An arbirary wavefront can be decomposed into a set of paraxial Gaussian beamlets (called "Gausslets" henceforth) and 
propagated through an optical system using raytracing. This method provides a convenient route to extend a geometric ray-tracing
algorithm to obtain a physical-optics modelling solution. The method has been successfully employed in a number of commercial
optical modelling packages, such as CODE-V, FRED and ASAP. Given a set of Gausslets, the vector (E) field can be evaluated at 
any position by summing the fields from each Gausslet. 

The method is accurate provided some restrictions are observed:

* The Gausslets must remain paraxial (i.e. more-or-less collimated).

* Individual Gausslets should interact with refracting or reflecting surfaces over a sufficiently small area such that the 
  surface may be considered locally parabolic.

Typically, if either of these restrictions are violated, the solution is to interrupt the ray-tracing operation to evaluate
the field, then decompose the field into a new set of Gausslets. In cases where the incident Gausslet is sampling too large
an area of the surface for it to be considered parabolic, the field is decomposed into a grid of smaller Gausslets. At the 
other extreme, Gausslets incident on a small aperture may be decomposed into a set of spatially overlapped Gausslets with
range of directions. Such an angle decomposition is sometimes known as a "Gabor Decomposition".

Raypier Gausslet Implementation
===============================

Raypier support Gausslet tracing through a GaussletCollection data-structure. Gausslets are an extension of the Ray data structure.
In fact, the dtype for a Gausslet object looks like::

    para_dtype = [('origin', '<f8', (3,)), ('direction', '<f8', (3,)), ('normal', '<f8', (3,)), ('length', '<f8')]

    gausslet_dtype = [('base_ray', ray_dtype), ('para_rays', para_dtype, (6,))]
    
I.e. we define a dtype for the parabasal rays which contains only the geometric information (origin, direction, length etc.) 
and omits any E-field polarisation or amplitude information. The gausslet is then composed of one `base_ray` (with 
regular ray_dtype) and 6 parabasal rays (with para_dtype).

Gausslets have their own set of predefined source-objects, found in :py:mod:`raypier.gausslet_sources`. 

Evaluating the E-field
======================

The nice thing about Gausslet ray-tracing is that you can evaluate the E-field at any point in your model. For script-based analysis,
you can give any GaussletCollection object (obtained from a source-object after a tracing operation) to the 
:py:func:`raypier.core.fields.eval_Efield_from_gausslets` function.
 
 .. py:module:: raypier.core.fields
 
Algorithms for evaluation of E-fields are provided here.
 
 .. py:function:: raypier.core.fields.eval_Efield_from_gausslets(gc : GaussletCollection, points : ndarray[N,3], wavelengths=None, blending=1.0) -> Efield ndarray[N,3]
                                                
    Calculates the vector E-field is each of the points given. The returned 
    array of field-vectors will have the same length as `points` and 
    has `numpy.complex128` dtype.
    
    :param GaussletCollection gc: The set of Gausslets for which the field should be calculated
    :param ndarray[N,3] points: An array of shape (N,3) giving the points at which the field will be evaluated.
    :param ndarray[] wavelengths: A 1d array containing the wavelengths to be used for the field calculation, 
                                    overriding the wavelengths data contained by the GaussletCollection object.
    :param float blending: The 1/width of each Gaussian mode at the evaluation points. A value of unity (the default),
                            means the parabasal rays are determined to be the 1/e point in the field amplitude.
                            
.. py:class:: EFieldSummation(gc : GaussletCollection, points : ndarray[N,3], wavelengths=None, blending=1.0) -> EFieldSummation object

    For situations where you wish to evaluate the E-field from a set of Gausslets with different sets of evaluation points,
    this class provides a small optimisation by performing the maths to convert ray-intercepts to Gaussian mode parameters
    up front.
    
    .. py:method:: evaluate(points : ndarray[N,3]) -> Efield ndarray[N,3]
    
        Called to calculate the E-field for the given points.
        
        
Beam Decomposition
==================

When a Beam-decomposition object intercepts a ray during the tracing operation, instead of immediately generating
child rays as most other `Traceable` objects do, the decomposition objects simply store the intercepted ray.
At the completion of tracing of the current ray generation (i.e. GaussletCollection), any decomposition-objects
which have received one or more rays then perform their decomposition-algorithm to generate a new set of rays to
be added to the other rays created in the last generation. The new rays created by the decomposition process
will, in general, not join originate from the end-point of the input-rays.


.. py:module:: raypier.gausslets

High level `Optics` objects for beam decomposition are provided here.

.. py:class:: raypier.gausslets.PositionDecompositionPlane(BaseDecompositionPlane)
    
    Defines a plane at which position-decomposition will be beformed.
    
    .. py:attribute:: radius
    
        Sets the radius used for capturing incoming rays. Rays outside of this will "miss"
        
    .. py:attribute:: curvature
    
        An approximate radius-of-curvature to the beam focus. This is used to improve the 
        phase-unwrapping of the wavefront. The default is zero, which means a plane-wave
        is assumed. Negative values imply a focus behind the decomposition plane (i.e.
        on the opposite side to the plane direction vector).
        
    .. py:attribute:: resolution
    
        Sets the resampling density of the decomposition, in terms of the number 
        of new rays per `radius` extent.
        
    .. py:attribute:: blending
    
        Sets the blending values for the new rays. The new rays will have Gaussian
        1/e**2 intensity widths equal to `spacing`/`blending`, where the `spacing`
        value is `radius`/`resolution`.
        
.. py:class:: AngleDecomposition(BaseDecompositionPlane)
 
    Defines a plane at which Gabor (angle)-decomposition is to be performed.
    
    .. py:attribute:: sample_spacing
    
        Sets the sample-spacing at the decomposition plane, in microns.
        
    .. py:attribute:: width
    
        A value in the range 1->512 to set the number of samples along the width of the sample-plane.
        
    .. py:attribute:: height
    
        A value in the range 1->512 to set the number of samples along the height of the sample-plane.
        
    .. py:attribute:: mask
    
        A 2d array with shape matching the (width, height) and dtype numpy.float64 . The array
        values should be in the range 0.0 -> 1.0. This will be used to mask the input E-field.
        
    .. py:attribute:: max_angle
    
        Limits the angular divergence of the outgoing rays.
        

.. py:module:: raypier.core.gausslets

The low-level beam-decomposition algorithms are found in this module. Two types of decomposition are available: position-decomposition
and angle-decomposition. Use the former when the Gausslets are found to be too wide at a particular surface in the optical path 
to re-sample the beam onto a set of more compact Gausslets. The later is used to simulate the effect of apertures much smaller than
the Gausslet widths, such that each Gausslet can be treated like a plane-wave and the field-distribution found using a 2d Fourier
transform.



    
    

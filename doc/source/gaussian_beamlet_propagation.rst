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

The nice thing about Gausslet ray-tracing is that you can evaluate the E-field at any point in your model. 

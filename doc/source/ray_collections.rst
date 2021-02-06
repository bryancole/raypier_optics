
=============================
RayCollection and Ray Objects
=============================

The most important object in Raypier is the RayCollection. This is (as you might guess from the name) a 1D array of Ray objects.

The Ray object (:py:class:`raypier.core.ctracer.Ray`) represents a single ray, and wraps an underlying C-structure. The single 
Ray object exists as a convenience for python scripting. The actual ray-tracing operation operates only RayCollection objects.

:py:class:`Ray` objects have the following attributes:

* origin - A 3-vector (tuple) of floats giving the start-point for the ray
* direction - A 3-vector giving the direction of the ray. This is always normalised to unit length
* E_vector - A 3-vector defining the polarisation-axes. This vector is unit-length and orthogonal to the `direction` vector
* normal - A 3-vector giving the unit-length surface normal of the face from which the ray originated. May be undefined for rays 
  which have not originated from a face intersection
* refractive_index - a complex value giving the refractive index of the material through which the ray has propagated
* E1_amp - the complex electric field amplitude for polarisations parallel to the E_vector axis
* E2_amp - the complex electric field amplitude for the polarisation orthogonal to the E_vector axis
* length - the geometric length of the ray from it's start-point to its termination at an intersecting face (or may be
  set to the `max length` of the ray, it no intersection has occured
* phase - An additional phase-factor that may be introduced  by face interactions. Currently, only used to hold the
  "grating phase" arising from diffraction-grating surfaces.
* accumulated_path - the total *optical* path length accumulated from the parent ray plus (parent length * real part of 
  refractive index)
* wavelength_idx - a index into the wavelength list (held by both the RayCollection and/or source object)
* parent_idx - the index of the parent-ray in the parent RayCollection object
* end_face_idx - the index into the Global Face List of the face at which the ray terminates (i.e. intersects). The Global Face
  List can be accessed on the `all_faces` atrtibure of the :py:class:`RayTraceModel` object.
* ray_type_idx - a bitfield indicating is the ray is a reflected or transmitted ray. Other ray-types may be defined in future

Rays have some additional read-only properties defined:

* power - the sum of squares of the E1_amp and E2_amp components
* amplitude - the square-root of the power
* jones_vector - Returns a 2-tuple (alpha, beta) representing the Jones Vector describing the polarisation state of the ray.
  See _https://en.wikipedia.org/wiki/Jones_calculus#Jones_vector
* E_left - Returns the complex electric field amplitude for the left circular polarisation state
* E_right - Returns the complex electric field amplitude for the right circular polarisation state
* ellipticity - Returns the ratio of the power in the right-hand polarisation state to the left-hand polarisation state
  I.e. A value of +1 indicates 100% right-circular polarisation, 0 indicates linear polarisation, -1 indicates 100% left
  polarisation.
* major_minor_axes - Returns a 2-tuple of unit-length vectors describing the major and minor polarisation axes

RayCollection objects have substantially the same attributes/properties as the Ray object, except that each property
returns a numpy array containing the values for all rays in the collection. 

RayCollection objects are iterable (yielding single Rays) and subscriptable. 

Creating RayCollections
-----------------------

If you need to create a RayCollection with a large number of rays (say, if you are writing your own Source class), the easiest 
method is to create a numpy array with the equivalent numpy dtype::

    >>> from raypier.api import ray_dtype, RayCollection
    >>> print(ray_dtype)
    [('origin', '<f8', (3,)), 
     ('direction', '<f8', (3,)), 
     ('normal', '<f8', (3,)), 
     ('E_vector', '<f8', (3,)), 
     ('refractive_index', '<c16'), 
     ('E1_amp', '<c16'), 
     ('E2_amp', '<c16'), 
     ('length', '<f8'), 
     ('phase', '<f8'), 
     ('accumulated_path', '<f8'), 
     ('wavelength_idx', '<u4'), 
     ('parent_idx', '<u4'), 
     ('end_face_idx', '<u4'), 
     ('ray_type_id', '<u4')]
     
Once you've created a numpy array with this dtype, you populate its fields as required. You can then create a RayCollection 
instance from this array using :py:meth:`RayCollection.from_array' classmethod. E.g.::

    >>> import numpy
    >>> my_rays = numpy.zeros(500, ray_dtype)
    >>> my_rays['direction'] = numpy.array([[0,1,1]],'d')
    <... assign other members as necessary ...>
    >>> rc = RayCollection.from_array(my_rays)

Note, data is always copied to and from RayCollections. The reason why we don't use memory views is that RayCollections have a
dynamic size and can grow in size by re-allocation of their memory. Numpy array, by contrast are static in size.

Likewise, we can convert a RayCollection to a numpy array using its :py:meth:`RayCollection.copy_as_array` method.::

    >>> arr = rc.copy_as_array()
    
    
   
    
    



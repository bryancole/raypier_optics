Cython branch: How to Add New Optics
====================================

When working with the cython raytrace branch, be sure you have the latest 
Cython version installed (as of writing, v0.12.1 was the latest).

Creating a new Face
-------------------

Subclasses of raytrace.ctracer.Face perform the mechanics of the ray-tracing operation.
Most existing subclasses are defined in cfaces.pyx. To define a new face object,
two "C" methods need to be defined:

    cdef double intersect_c(vector_t p1, vector_t p2)

Computes the intersection of a ray with the surface, where p1 and p2 are two 
points defining the start and end points of an incoming
ray, in the local coordinate system. The method returns a single floating-point 
value representing the distance along the ray where an intersection occurs. If
no intersection is found, a value <= 0 may be returned (I tend to return zero,
if no intersection occurs).

    cdef vector_t compute_normal_c(vector_t p)

Computes the outward normal vector for the surface at the given point, p. These
are also defined w.r.t. the local coordinate system.

Note, these two methods are 'cdef-methods' and hence not callable directly from 
python. Be sure to type all the variables you use, to be sure of good performance.
Face objects also have python-callable methods compute_normal() and intersect(). 
These call the cdef methods internally. Don't bother trying to overload these
from python as it won't work (and even if it did, performance would be horrible).

The vector_t structure is used to represent 3D points (i.e. it's a 3d vector),
and has .x and .y and .z members. Outside the comfort of python/numpy, vector
algebra is significantly less pleasant. The ctracer module defines a number of
inline functions to (partially) simplify working with vector_t's. For vector 
arithmatic, you can use:
 - addvv_(vector_t a, vector_t b) #add two vectors
 - addvs_(vector_t a, double b) #adds a vector and a scalar
 - multvv_(vector_t a, vector_t b) #multiplies two vectors
 - multvs_(vector_t a, double b) #multiplies a vector and a scalar
 etc. 
 Similar function are defined for subtraction and element-wise division. The
'vv' function denote operations on two vectors, 'vs' denotes operations on a
vector and a scalar. The trailing underscore is a convention I'm adopting to
indicate that the method is a 'C'-method (i.e. not python callable).

Some of the other utility functions are:

  cdef inline vector_t set_v(vector_t v, object O)
  
    Takes a vector v and an 3-sequence O. Copies the values from the sequence
    into the vector and returns it (why does it need the v argument, though?).

  cdef inline double sep_(vector_t p1, vector_t p2)

    Computed the linear separation between two points

  cdef inline double dotprod_(vector_t a, vector_t b)

    Calculates the dot-product of the given vectors

  cdef inline vector_t cross_(vector_t a, vector_t b)

    Calculates the vector-product of it's arguments

  cdef inline vector_t norm_(vector_t a)

    Computes the normalised vector from its input (i.e. scales it's input
    to unit magnitude).

  cdef inline double mag_(vector_t a)
  
    Calculates the magnitude of it's input
  
  cdef inline double mag_sq_(vector_t a)

    Calculates the square of the magnitude of it's input

  cdef inline vector_t invert_(vector_t v)

    Inverts it's input


Most of the cdef functions and cdef-methods have companion def functions/methods
which are callable from python and wrap the associated cdef operation. These
are mostly used in testing. The python-callable versions incure the normal python
interpreter overhead, and hence are not directly useful in the fast C-level 
tracing process.

The intersect_c and compute_normal_c cdef methods are all you need to make a 
new face. Obviously, you can add additional methods / attributes to implement
the functionality required.

You can optionally give your custom Face class a params attribute (define it as 
a class-attribute) which should contain a list of attribute names which should
be synchronised from the face owner to the face when a tracing operation is 
initiated.

Creating a new Traceable
------------------------

Traceables (i.e. subclasses of raytrace.bases.Traceable) are the basic unit of
an optical model. Most of the functionality in the Traceable subclasses is 
to implement their VTK visual representation. The ray-tracing operation is 
handled by Face objects. Traceables own an instance of a ctracer.FaceList 
which turn has a list of Face objects (I choose FaceList *has a* list of faces,
rather than FaceList *is a* list of faces, as it was simpler to implement). 
The FaceList contains the coordinate transform which maps betwee global coords 
and the local coords of the Traceable. Thus, all Faces belonging to a Traceable
share a common transform.

To create a new Traceable, you subclass Traceable or some other more suitable
subclass (transmitting optical components can derive from raytrace.bases.Optic,
which provides a complex refractive index). You should define a new _faces_default
method which creates the FaceList for that object and populates it with the
Faces appropriate to the object. Simple synchronisation between the Traceable
and the Faces can be handled using the params Face class attribute described
above. In most cases, more sophisticated synchronisation is required and can
be handled using trait-notifications for all traits on which the Faces depends.

The physics of ray-scattering (i.e. the generation of new rays at the point of
intersection) is handled by ctracer.InterfaceMaterial objects. InterfaceMaterial
is an abstract base class. There are two concrete subclasses defined in the ctracer
module: PECMaterial and DielectricMaterial. The former represents a perfect metal
reflector. The later is a normal dielectric surface. Typically, an Optic (or
Traceable subclass) will have an InterfaceMaterial trait. This will be passed
to it's faces in the _faces_default method (so all faces share the same
InterfaceMaterial). However, this is not a requirement: faces can have independent 
materials, or share them.

Custom Interface Materials
--------------------------

InterfaceMaterial subclasses provide a cdef method 

  cdef eval_child_ray_c(self, ray_t *old_ray, 
                            unsigned int ray_idx, 
                            vector_t point, vector_t normal,
                            RayCollection new_rays)
                            
This is called for each ray intersection to create a new ray. The arguments
are as follows:

  old_Ray - a pointer to the incoming ray_t structure
  ray_idx - the index of the incoming ray in it's RayCollection array
  point - the position, in global coords, of the intersection
  normal - the normal vector of the surface, at the point of intersection
  new_rays - the target RayCollection for new rays
  
This method should call new_rays.add_new_ray() to create as many new rays as
necessary. Thus, multiple ray generation can occur at an intersection (as might
be found for a diffracting interface material).

Cython Tips and Tricks
----------------------

If you find performance is less than you expected, try running "cython -a yourfile.pyx"
(replace yourfile.pyx with whatever .pyx file you're analysing, cfaces.pyx maybe).
This produces a html-version of your file, with highlighting to show where the python 
API is being invoced. The less yellow the better (and red-highlights indicate
real performance bottlenecks). This is a *very* cool feature of Cython.

Avoid cpdefs (i.e. methods with automatically created python wrappers), as extra
overhead is incured to check for python overloading.

Creating and destroying python objects is expensive (it invokes the garbage
collector / changes ref-counts etc.). However, read-only access to items in 
lists is fast.

Surprisingly, I can find no speed benefit in passing parameters by reference, 
compared to passing by values (for fixed-size types, at least).

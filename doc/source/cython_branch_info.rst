New Cython-base raytracing framework
====================================

The core tracing functions of Raytrace have been re-written using Cython, a 
tool for creating compiled C-extension modules using a python-like language. The
new architecture performs substantially better than the previous numpy-based design
(In a 4-mirror model, 8x faster than numpy for 60 rays, 25x faster for 1000 rays
and 75x faster for 2500 rays).

The main tracing functionality is provided by three new (compiled) modules: 
 - ctracer
 - cfaces
 - cmaterials
 
"ctracer" provides the following extension types:
 - Ray
 - RayCollection
 - Face (abstract base class)
 - FaceList
 - InterfaceMaterial (abstract base class)
 - PECMaterial
 - DielectricMaterial
 - Transform
 
 
.. module:: ctracer
..moduleauthor:: Bryan Cole <bryancole.cam@googlemail.com>

 
.. class:: Ray

  Ray is a wrapper class around the ray_t structure.
  
.. class:: RayCollection(max_size)

  A RayCollection is a wrapper around a C-array of ray_t.
  
  .. method:: add_ray(ray)
  
    ray is an instance of Ray. Adds the given Ray to the array
  
  
.. module:: cfaces

  This module contains the concrete implementation of Face classes
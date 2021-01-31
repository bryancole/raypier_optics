
from .ctracer import trace_segment, trace_gausslet, RayCollection

from itertools import chain
import numpy


def trace_rays(input_rays, face_lists, recursion_limit=100, max_length=100.0):
    """
    Core ray-tracing routine. Takes either a RayCollection or GaussletCollection
    and traces the rays non-sequentially through the given list of FaceList objects.
    
    The input_rays should already have a consistent wavelengths property set.
    
    returns - (traced_rays, all_faces)
            where traced_rays is a list of RayCollection or GaussletCollections 
            representing the sequence of ray generations. The 'all_faces' list
            is a list of cfaces.Face objects. The end_face_idx member of each
            traced ray indexes into this list to give the face where it terminates.
    """
    input_rays.reset_length(max_length)
    traced_rays = []
    trace_func = trace_segment if isinstance(input_rays, RayCollection) else trace_gausslet
    count = 0
    wavelengths = numpy.asarray(input_rays.wavelengths)
    
    all_faces = list(chain(*(fs.faces for fs in face_lists)))
    for i, f in enumerate(all_faces):
        f.idx = i
        f.count = 0 #reset intersection count
        f.update() #This was a bad idea.
        f.material.wavelengths = wavelengths
        f.max_length = max_length
    
    face_sets = list(face_lists)
    decomp_faces = [f for f in all_faces if f.material.is_decomp_material()]
    
    rays = input_rays
    while rays.n_rays>0 and count<recursion_limit:
        traced_rays.append(rays)
        rays = trace_func(rays, face_sets, all_faces, 
                                     max_length=max_length,
                                     decomp_faces=decomp_faces)
        count += 1
    
    return traced_rays, all_faces 


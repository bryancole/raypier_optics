#!/usr/bin/env python
from cfaces cimport FaceList, Intersection

cdef extern from "math.h":
    double sqrt

cdef double sep(double[3] p1, double[3] p2):
    double val
    val = 0
    for i in range(3):
        val += (p2[i] - p1[i])**2
    return sqrt(val)

cdef object trace_segment(object rays, object optics):
    FaceList face_set
    int size, i
    float[3] P1, P2
    float max_length
    Intersection inter, inter2
    
    face_sets = [o.faces for o in optics]
    size = rays.number
    origin = rays.origin
    direction = rays.direction
    max_length = rays.max_length
    
    #need to allocate the output rays here
    
    for i in range(size):
        P1 = origin[i]
        P2 = P1 + (direction[i] * max_length)
        inter = face_sets[0].intersect(P1, P2)
        dist = sep(inter.point, P1)
        for face_set in face_sets[1:]: #need to make this a C array
            inter2 = face_set.intersect(P1, P2)
            dist2 = sep(inter2.point, P1)
            if dist2 < dist:
                dist = dist2
                inter = inter2
        
        #now compute normal
        normal = inter.face.compute_normal(inter.point)
        
        #evaluate new ray
        
                
        
        
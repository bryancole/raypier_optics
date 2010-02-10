#!/usr/bin/env python

cdef extern from "math.h":
    double sqrt(double arg)
    double INFINITY
    
from stdlib cimport malloc, free

############################################
### C type declarations for internal use ###
############################################

cdef struct vector_t:
    double x,y,z
    
cdef struct complex_t:
    double real
    double imag
    
cdef struct ray_t:
    #vectors
    vector_t origin, direction, normals, E_vector
    #complex attribs
    complex_t refractive_index, E1_amp, E2_amp
    #simple attribs
    double length, wavelength
    #reference ids to related objects
    unsigned int parent_idx, end_face_idx
    ##objects
    #object face, end_face, child_refl, child_trans
    
cdef struct transform_t:
    double m00, m01, m02, m10, m11, m12, m20, m21, m22
    double tx, ty, tz
    
    
cdef struct intersection_t:
    unsigned int face_idx #the intersecting face
    vector_t point #position of intersection
    double dist #fractional of ray between origin and intersection


##############################
### Vector maths functions ###
##############################

cdef inline vector_t transform_c(transform_t t, vector_t p)
cdef inline vector_t set_v(vector_t v, object O)
cdef inline double sep_(vector_t p1, vector_t p2)
cdef inline vector_t multvv_(vector_t a, vector_t b)
cdef inline vector_t multvs_(vector_t a, double b)
cdef inline vector_t addvv_(vector_t a, vector_t b)
cdef inline vector_t addvs_(vector_t a, double b)
cdef inline vector_t subvv_(vector_t a, vector_t b)
cdef inline vector_t subvs_(vector_t a, double b)

##################################
### Python extension types
##################################

cdef class Transform:
    cdef transform_t trans

cdef class Intersection:
    cdef intersection_t inter

cdef class Ray:
    cdef ray_t ray

cdef class RayCollection:
    cdef ray_t *rays
    cdef readonly unsigned long n_rays, max_size
    cdef public RayCollection parent
        
    cdef add_ray_c(self, ray_t r)
    
    
cdef class Face(object):
    cdef public object owner
    cdef public char *name
    cdef public double tolerance
    cdef public int idx #index in the global face list
    
    cdef vector_t intersect_c(self, vector_t p1, vector_t p2)

    cdef vector_t compute_normal_c(self, vector_t p)
    
    cdef ray_t eval_child_ray_c(self, ray_t old_ray, int ray_idx, 
                                vector_t p)


cdef class FaceList(object):
    """A group of faces which share a transform"""
    cdef transform_t trans
    cdef transform_t inv_trans
    cdef public list faces
     
    cdef intersection_t intersect_c(self, vector_t P1, vector_t P2, double max_length)

##################################
### Python module functions
##################################

cdef RayCollection trace_segment_c(RayCollection rays, 
                                    list face_sets, 
                                    list all_faces)

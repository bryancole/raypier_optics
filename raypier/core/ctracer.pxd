#!/usr/bin/env python
from _operator import invert

cdef extern from "float.h":
    #double INFINITY
    double DBL_MAX

cdef:
    double INF
    int NPARA=6


from libc.stdlib cimport malloc, free
cimport numpy as np_

############################################
### C type declarations for internal use ###
############################################

cdef struct vector_t:
    double x,y,z

cdef struct orientation_t:
    vector_t normal, tangent


IF UNAME_SYSNAME == "Windows":
    ctypedef double complex complex_t
ELSE:
    ctypedef double complex complex_t
#     cdef struct complex_t:
#         double real
#         double imag
        
cdef:
    public unsigned int REFL_RAY=1<<0
    public unsigned int GAUSSLET=1<<1
    public unsigned int PARABASAL=1<<2


cdef packed struct ray_t:
    #vectors
    vector_t origin, direction, normal, E_vector
    #complex attribs
    complex_t refractive_index, E1_amp, E2_amp
    #simple attribs
    double length, phase, accumulated_path
    #reference ids to related objects
    unsigned int wavelength_idx, parent_idx, end_face_idx, ray_ident, ray_type_id
    
    ##objects
    #object face, end_face, child_refl, child_trans
    
### Structure for parabasal rays. These only carry geometric information.
cdef packed struct para_t:
    #vectors
    vector_t origin, direction, normal
    #simple attribs
    double length

    
### A Gausslet combines a base-ray with 6 parabasal rays
cdef packed struct gausslet_t:
    ray_t base_ray
    para_t para[6]
    

cdef struct transform_t:
    double m00, m01, m02, m10, m11, m12, m20, m21, m22
    double tx, ty, tz

cdef struct ray_pair_t:
    ray_t trans, refln
    
cdef struct intersect_t:
    double dist
    int face_idx
    int piece_idx
    

##############################
### Vector maths functions ###
##############################

cdef  vector_t transform_c(transform_t t, vector_t p) nogil
cdef  vector_t rotate_c(transform_t t, vector_t p) nogil
cdef  vector_t set_v(object O)
cdef  double sep_(vector_t p1, vector_t p2) nogil
cdef  vector_t multvv_(vector_t a, vector_t b) nogil
cdef  vector_t multvs_(vector_t a, double b) nogil
cdef  vector_t addvv_(vector_t a, vector_t b) nogil
cdef  vector_t addvs_(vector_t a, double b) nogil
cdef  vector_t subvv_(vector_t a, vector_t b) nogil
cdef  vector_t subvs_(vector_t a, double b) nogil
cdef  double dotprod_(vector_t a, vector_t b) nogil
cdef  vector_t cross_(vector_t a, vector_t b) nogil
cdef  vector_t norm_(vector_t a) nogil
cdef  double mag_(vector_t a) nogil
cdef  double mag_sq_(vector_t a) nogil
cdef  vector_t invert_(vector_t v) nogil


##################################
### Python extension types
##################################

cdef class Transform:
    cdef transform_t trans
    

cdef class Ray:
    cdef:
        ray_t ray
        double max_length
    

cdef class ParabasalRay:
    cdef:
        para_t ray
        double max_length
    
    
cdef class Gausslet:
    cdef gausslet_t gausslet
    
    
cdef class RayArrayView:
    cdef void set_ray_c(self, unsigned long i, ray_t ray) nogil
    cdef ray_t get_ray_c(self, unsigned long i) nogil
    cdef unsigned long get_n_rays(self) nogil
    
    
cdef class RayAggregator:
    cdef:
        ray_t **rays
        readonly unsigned long n_threads
        unsigned long[:] n_rays
        unsigned long[:] max_size
        
    cdef void add_ray_c(self, ray_t r) nogil
    cdef clear_c(self)
    cdef get_aggregated_c(self)
    

cdef class RayCollection(RayArrayView):
    cdef: 
        ray_t *rays
        readonly unsigned long n_rays, max_size
        RayCollection _parent
        double[:] _wavelengths
        
        int[:,:] _neighbours
        double _mtime        

    cdef add_ray_c(self, ray_t r) nogil
    cdef void reset_length_c(self, double max_length)
    
    cdef double get_mtime(self, unsigned long guard)
    
    cdef void _eval_neighbours(self, int[:,:] pnb)
    
    
cdef class GaussletCollection:
    cdef: 
        gausslet_t *rays
        readonly unsigned long n_rays, max_size
        GaussletCollection _parent
        double[:] _wavelengths

    cdef add_gausslet_c(self, gausslet_t r)
    cdef void extend_c(self, GaussletCollection gc)
    cdef void reset_length_c(self, double max_length)
    
    
cdef class GaussletBaseRayView(RayArrayView):
    cdef:
        GaussletCollection owner
    

cdef class RayCollectionIterator:
    cdef:
        RayCollection rays
        unsigned int counter
        
        
cdef class GaussletCollectionIterator:
    cdef:
        GaussletCollection rays
        unsigned int counter


cdef class InterfaceMaterial(object):
    """Abstract base class for objects describing
    the materials characterics of a Face
    """
    cdef double[:] _wavelengths

    cdef void eval_child_ray_c(self, ray_t *old_ray,
                            unsigned int ray_idx,
                            vector_t point,
                            orientation_t orient,
                            void* new_rays) nogil
    
    cdef para_t eval_parabasal_ray_c(self, ray_t *base_ray,  
                                     vector_t direction, #incoming ray direction
                                   vector_t point, #position of intercept
                                   orientation_t orient,
                                   unsigned int ray_type_id, #indicates if it's a transmitted or reflected ray 
                                   ) nogil
    
    cdef void eval_decomposed_rays_c(self, GaussletCollection child_rays)

    cdef on_set_wavelengths(self)


cdef class Distortion:
    cdef vector_t z_offset_and_gradient_c(self, double x, double y) nogil
    cdef double z_offset_c(self, double x, double y) nogil
    
    
cdef class Shape:
    cdef bint point_inside_c(self, double x, double y) nogil


cdef class Face(object):
    cdef:
        public object owner
        public char *name
        public double tolerance
        public int idx #index in the global face list
        public double max_length
        public InterfaceMaterial material
        public short int invert_normal
        public unsigned int count    

    cdef intersect_t intersect_c(self, vector_t p1, vector_t p2, int is_base_ray) nogil

    cdef vector_t compute_normal_c(self, vector_t p, int piece) nogil
    cdef vector_t compute_tangent_c(self, vector_t p, int piece) nogil


cdef class FaceList(object):
    """A group of faces which share a transform"""
    cdef:
        transform_t trans
        transform_t inv_trans
        public np_.ndarray _faces
        void** _faces_ptr
        public object owner
        int size

    cpdef void sync_transforms(self)
    cdef intersect_t intersect_c(self, ray_t *ray, vector_t end_point) nogil
    cdef intersect_t intersect_one_face_c(self, ray_t *ray, vector_t end_point, int face_idx) nogil
    cdef int intersect_para_c(self, para_t *ray, vector_t ray_end, Face face)
    cdef orientation_t compute_orientation_c(self, Face face, vector_t point, int piece) nogil


##################################
### Python module functions
##################################

cdef RayCollection trace_segment_c(RayCollection rays,
                                    np_.ndarray face_sets,
                                    np_.ndarray all_faces,
                                    list decomp_faces,
                                    float max_length)

cdef RayCollection trace_one_face_segment_c(RayCollection rays,
                                    FaceList face_set,
                                    int face_idx,
                                    list all_faces,
                                    list decomp_faces,
                                    float max_length)

cdef GaussletCollection trace_gausslet_c(GaussletCollection gausslets, 
                                    list face_sets, 
                                    list all_faces,
                                    list decomp_faces,
                                    double max_length)

cdef GaussletCollection trace_one_face_gausslet_c(GaussletCollection gausslets, 
                                    FaceList face_set,
                                    int face_idx, 
                                    list all_faces,
                                    list decomp_faces,
                                    double max_length)

cdef double ray_power_(ray_t ray)

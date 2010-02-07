#!/usr/bin/env python

#cdef extern from "math.h":
#    double sqrt
#    double INFINITY
    
from stdlib cimport malloc, free

    
ctypedef double vector[3]

#cdef double sep(vector p1, vector p2):
#    cdef double val
#    val = 0
#    for i in range(3):
#        val += (p2[i] - p1[i])**2
#    return sqrt(val)
    
cdef struct complex_t:
    double real
    double imag
    
cdef struct ray_t:
    #vectors
    vector origin, direction, normals, E_vector
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
    vector point #position of intersection
    double dist #fractional of ray between origin and intersection


cdef class Ray:
    cdef ray_t ray
    
    def __cinit__(self, **kwds):
        for k in kwds:
            setattr(self, k, kwds[k])
            
    def __repr__(self):
        return "Ray(o=%s, d=%s)"%(str(self.origin),
                                            str(self.direction))
                
    property origin:
        """Origin coordinates of the ray"""
        def __get__(self):
            return (self.ray.origin[0],self.ray.origin[1],self.ray.origin[2])
        
        def __set__(self, v):
            self.ray.origin[0] = v[0]
            self.ray.origin[1] = v[1]
            self.ray.origin[2] = v[2]
            
    property direction:
        """direction of the ray, normalised to a unit vector"""
        def __get__(self):
            return (self.ray.direction[0],self.ray.direction[1],self.ray.direction[2])
        
        def __set__(self, v):
            self.ray.direction[0] = v[0]
            self.ray.direction[1] = v[1]
            self.ray.direction[2] = v[2]
            
    property normals:
        """normal vector for the face which created this ray"""
        def __get__(self):
            return (self.ray.normals[0],self.ray.normals[1],self.ray.normals[2])
        
        def __set__(self, v):
            self.ray.normals[0] = v[0]
            self.ray.normals[1] = v[1]
            self.ray.normals[2] = v[2]
            
    property E_vector:
        """Unit vector, perpendicular to the ray direction,
        which gives the direction of E-field polarisation"""
        def __get__(self):
            return (self.ray.E_vector[0],self.ray.E_vector[1],self.ray.E_vector[2])
        
        def __set__(self, v):
            self.ray.E_vector[0] = v[0]
            self.ray.E_vector[1] = v[1]
            self.ray.E_vector[2] = v[2]
            
    property length:
        """The length of the ray. This is infinite in 
        unterminated rays"""
        def __get__(self):
            return self.ray.length
        
        def __set__(self, double v):
            self.ray.length = v
        
    property refractive_index:
        """complex refractive index through which
        this ray is propagating"""
        def __get__(self):
            return complex(self.ray.refractive_index.real,
                            self.ray.refractive_index.imag)
        def __set__(self, v):
            self.ray.refractive_index.real = v.real
            self.ray.refractive_index.imag = v.imag
            
    property E1_amp:
        """Complex amplitude of the electric field polarised
        parallel to the E_vection."""
        def __get__(self):
            return complex(self.ray.E1_amp.real,
                            self.ray.E1_amp.imag)
        def __set__(self, v):
            self.ray.E1_amp.real = v.real
            self.ray.E1_amp.imag = v.imag
            
    property E2_amp:
        """Complex amplitude of the electric field polarised
        perpendicular to the E_vection"""
        def __get__(self):
            return complex(self.ray.E2_amp.real,
                            self.ray.E2_amp.imag)
        def __set__(self, v):
            self.ray.E2_amp.real = v.real
            self.ray.E2_amp.imag = v.imag
            
    property parent_idx:
        """Index of the parent ray in parent RayCollection
        """
        def __get__(self):
            return self.ray.parent_idx
        
        def __set__(self, int v):
            self.ray.parent_idx = v
            
    property end_face_idx:
        """Index of the terminating face, in the global
        face list (created for each tracing operation)
        """
        def __get__(self):
            return self.ray.end_face_idx
        
        def __set__(self, int v):
            self.ray.end_face_idx = v


cdef class RayCollection:
    cdef ray_t *rays
    cdef readonly unsigned long n_rays, max_size
    cdef public RayCollection parent
    
    def __cinit__(self, size_t max_size):
        self.rays = <ray*>malloc(max_size*sizeof(ray_t))
        self.n_rays = 0
        self.max_size = max_size
        
    def __dealloc__(self):
        free(self.rays)
        
    cdef add_ray_c(self, ray_t r):
        self.rays[self.n_rays] = r
        self.n_rays += 1
        
    def add_ray(self, Ray r):
        if self.n_rays == self.max_size-1:
            raise ValueError("RayCollection is full up")
        self.rays[self.n_rays] = r.ray
        self.n_rays += 1
        
    def add_ray_list(self, list rays):
        cdef int i
        if self.n_rays + len(rays) >= self.max_size:
            raise IndexError("RayCollection size is too small to hold this many rays")
        for i in xrange(len(rays)):
            if not isinstance(rays[i], Ray):
                raise TypeError("ray list contains non-Ray instance at index %d"%i)
        for i in xrange(len(rays)):
            self.add_ray_c((<Ray>rays[i]).ray)
        
    def clear_ray_list(self):
        self.n_rays = 0
        
    def get_ray_list(self):
        cdef int i
        cdef list ray_list = []
        cdef Ray r
        for i in xrange(self.n_rays):
            r = Ray()
            r.ray = self.rays[i]
            ray_list.append(r)
        return ray_list
    
    def __getitem__(self, int idx):
        cdef Ray r
        if idx >= self.n_rays:
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        r = Ray()
        r.ray = self.rays[idx]
        return r
    
    def __setitem__(self, int idx, Ray r):
        if idx >= self.n_rays:
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        self.rays[idx] = r.ray
    
    
cdef class Face(object):
    cdef public object owner
    cdef public char *name
    cdef public double tolerance
    cdef public int idx #index in the global face list
    
    def __cinit__(self):
        self.name = "base Face class"
        self.tolerance = 0.0001
    
    cdef double intersect_c(self, vector p1, vector p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        return 0.5

#    cdef vector compute_normal_c(self, vector p):
#        return Vector(p.z,p.y,p.x)
    
    cdef ray_t eval_child_ray_c(self, ray_t old_ray, int ray_idx, 
                                vector p):
        return ray


cdef class FaceList(object):
    """A group of faces which share a transform"""
    cdef transform_t transform
    cdef transform_t inverse_transform
    cdef public list faces
     
    cdef intersection_t intersect_c(self, vector P1, vector P2, double max_length):
        """Finds the face with the nearest intersection
        point, for the ray defined by the two input points,
        P1 and P2 (in global coords).
        """
        cdef vector p1, p2, point
        cdef list faces
        cdef double d, dist=INFINITY
        cdef Face nearest=None
        cdef Intersection inter
        
        p1 = transform(self.transform, P1)
        p2 = transform(self.transform, P2)
        
        faces = self.faces
        
        for i in xrange(len(faces)):
            point = (<Face>(faces[i])).intersect_c(p1, p2)
            d = sep(p1, point)
            if 0.0 < d < dist:
                dist = d
                nearest = <Face>(faces[i])
        
        inter.face = nearest
        inter.point = transform(self.inverse_transform, point)
        inter.dist = dist/sep(p1,p2)
        return inter


cdef RayCollection trace_segment_c(RayCollection rays, 
                                    list face_sets, 
                                    list all_faces):
    cdef:
        FaceList face_set #a FaceList
        int size, i
        vector P1, P2, normal
        float max_length
        intersection_t inter, inter2
        int n_sets=len(face_sets)
        ray_t ray
   
    #need to allocate the output rays here
    new_rays = RayCollection(rays.n_rays)
   
    for i in range(size):
        ray = rays.rays[i]
        P1 = ray.origin
        P2 = P1 + (ray.direction * max_length)
        inter = (<FaceList>(face_sets[0])).intersect_c(P1, P2)
        dist = sep(inter.point, P1)
        for j in xrange(n_optics-1):
            face_set = face_sets[j+1]
            inter2 = (<FaceList>face_set).intersect_c(P1, P2)
            dist2 = sep(inter2.point, P1)
            if dist2 < dist:
                dist = dist2
                inter = inter2
                
        if dist <= INFINITY:
            face = all_faces[inter.face_idx]
            #evaluate new ray
            new_rays.add_ray_c(face.eval_child_ray_c(ray, i, 
                                                    inter.point))
    return new_rays

        
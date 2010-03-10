#!/usr/bin/env python

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

"""
Cython module for Face definitions
"""
cdef extern from "math.h":
    double INFINITY
    double sqrt(double)

from ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_


cdef class CircularFace(Face):
    cdef public double diameter, offset, z_plane
    
    params = ['diameter', 'offset']
    
    def __cinit__(self, **kwds):
        self.z_plane = kwds.get('z_plane', 0.0)
    
    cdef int intersect_c(self, vector_t p1, vector_t p2, ray_t *ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          ray - pointer to the input ray
          
        returns:
          idx - -1 if no intersection is found *or* the distance to the 
                intersection is larger than the existing ray.length. OTherwise,
                this is set to the intersecting face idx
        """
        cdef:
            double max_length, h, X, Y, d=self.diameter
            
        #print "CFACE", p1, p2
        max_length = sep_(p1, p2)
        h = (self.z_plane-p1.z)/(p2.z-p1.z)
        if (h<self.tolerance) or (h>1.0):
            #print "H", h
            return -1
        X = p1.x + h*(p2.x-p1.x) - self.offset
        Y = p1.y + h*(p2.y-p1.y)
        if (X*X + Y*Y) > (d*d/4):
            #print "X", X, "Y", Y
            return -1
        h *= max_length
        if h > ray.length:
            return -1
        ray.length = h
        ray.end_face_idx = self.idx
        return self.idx

    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef vector_t normal
        normal.x=0
        normal.y=0
        normal.z=-1
        return normal
    
    
cdef class SphericalFace(Face):
    cdef public double diameter, curvature, z_height
    
    params = ['diameter', 'curvature']
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_plane', 0.0)
    
    cdef int intersect_c(self, vector_t r, vector_t p2, ray_t *ray):
        """Intersects the given ray with this face.
        
        params:
          r - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          ray - pointer to the input ray
          
        returns:
          idx - -1 if no intersection is found *or* the distance to the 
                intersection is larger than the existing ray.length. OTherwise,
                this is set to the intersecting face idx
        """
        cdef:
            double A,B,C,D, cz
            vector_t s, d, pt
            
        s = subvv_(p2, r)
        cz = self.z_height - self.curvature
        d = r
        d.z -= cz
        
        A = mag_sq_(s)
        B = 2*dotprod_(s,d)
        C = mag_sq_(d) - self.curvature**2
        D = B*B - 4*A*C
        
        if D < 0: #no intersection with sphere
            return -1
        
        D = sqrt(D)
        
        A = (-B+D)/(2*A) #1st root
        pt = addvv_(r, multvs_(s, A))
        if pt.z < cz: #wrong root
            A = (-B-D)/(2*A) #2nd root
            pt = addvv_(r, multvs_(s, A))
        
        if A>1.0 or A<self.tolerance:
            return -1
        if (pt.x*pt.x + pt.y*pt.y) > (self.diameter*self.diameter/4.):
            return -1
        
        A *= sep_(r, p2)
        if A > ray.length:
            return -1
        ray.length = A
        ray.end_face_idx = self.idx
        return self.idx
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        p.z -= (self.z_height - self.curvature)
        return norm_(p)
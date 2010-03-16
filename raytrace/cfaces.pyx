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

import numpy as np
cimport numpy as np


cdef class CircularFace(Face):
    cdef public double diameter, offset, z_plane
    
    params = ['diameter', 'offset']
    
    def __cinit__(self, **kwds):
        self.z_plane = kwds.get('z_plane', 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2):
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
            double max_length = sep_(p1, p2)
            double h = (self.z_plane-p1.z)/(p2.z-p1.z)
            double X, Y, d=self.diameter
            
        #print "CFACE", p1, p2
        
        if (h<self.tolerance) or (h>1.0):
            #print "H", h
            return 0
        X = p1.x + h*(p2.x-p1.x) - self.offset
        Y = p1.y + h*(p2.y-p1.y)
        if (X*X + Y*Y) > (d*d/4):
            #print "X", X, "Y", Y
            return 0
        return h * max_length

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
    
    cdef double intersect_c(self, vector_t r, vector_t p2):
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
            return 0
        
        D = sqrt(D)
        
        A = (-B+D)/(2*A) #1st root
        pt = addvv_(r, multvs_(s, A))
        if pt.z < cz: #wrong root
            A = (-B-D)/(2*A) #2nd root
            pt = addvv_(r, multvs_(s, A))
        
        if A>1.0 or A<self.tolerance:
            return 0
        if (pt.x*pt.x + pt.y*pt.y) > (self.diameter*self.diameter/4.):
            return 0
        return A * sep_(r, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        p.z -= (self.z_height - self.curvature)
        return norm_(p)
    
    
cdef class ExtrudedPlanarFace(Face):
    cdef:
        double z1, z2, x1_, y1_, x2_, y2_
        vector_t normal
        
    property x1:
        def __get__(self):
            return self.x1_
        
        def __set__(self, double v):
            self.x1_ = v
            self.calc_normal()
            
    property y1:
        def __get__(self):
            return self.y1_
        
        def __set__(self, double v):
            self.y1_ = v
            self.calc_normal()
            
    property x2:
        def __get__(self):
            return self.x2_
        
        def __set__(self, double v):
            self.x2_ = v
            self.calc_normal()
            
    property y2:
        def __get__(self):
            return self.y2_
        
        def __set__(self, double v):
            self.y2_ = v
            self.calc_normal()
            
    cdef calc_normal(self):
        cdef vector_t n
            
        n.y = self.x2 - self.x1
        n.x = self.y1 - self.y2
        n.z = 0
        self.normal = norm_(n)
        
    cdef double intersect_c(self, vector_t r, vector_t p2):
        cdef: 
            vector_t s, u, v
            double a, dz
            
        u.x = self.x1
        u.y = self.y1
        
        v.x = self.x2 - u.x
        v.y = self.y2 - u.y
        
        s = subvv_(p2, r)
        
        #fractional distance of intersection along edge
        a = (s.y*(u.x-r.x) - s.x*(u.y-r.y)) / (s.x*v.y - s.y*v.x)
        if a<0:
            return 0
        if a>1:
            return 0
        #distance of intersection along ray (in XY plane)
        a = (v.x*(r.y-u.y) - v.y*(r.x-u.x)) / (s.x*v.y - s.y*v.x)
        
        #distance in 3D
        dz = a*(p2.z - r.z)
        
        if self.z1 < (r.z+dz) < self.z2:
            return sqrt(a*a + dz*dz)
        else:
            return 0
    
    cdef vector_t compute_normal_c(self, vector_t p):
        return self.normal
    
    
cdef point_in_polygon_c(double X, double Y,  obj):
    cdef int i, size, ct
    cdef double y1, y2, h, x, x1, x2
    cdef np.ndarray[np.float64_t, ndim=2] pts
    
    pts = obj
    size = pts.shape[0]
    
    y1 = pts[size-1,1]
    x1 = pts[size-1,0]
    for i in xrange(size):
        y2 = pts[i,1]
        x2 = pts[i,0]
        h = (Y - y1) / (y2 - y1)
        if 0 < h <= 1.0:
            x = x1 + h*(x2 - x1)
            if x > X:
                ct += 1
        y1 = y2
        x1 = x2
    if ct%2:
        return 0
    else:
        return 1
    
    
def point_in_polygon(double X, double Y, point_list):
    pts = np.ascontiguousarray(point_list, dtype=np.float64)
    assert pts.shape[1]==2
    return bool(point_in_polygon_c(X, Y, pts))
    
    
cdef class PolygonFace(Face):
    cdef public double z_plane
    cdef object _xy_points
    
    property xy_points:
        def __get__(self):
            return self._xy_points
        
        def __set__(self, pts):
            data = np.ascontiguousarray(pts, dtype=np.float64)
            assert data.shape[1]==2
            self._xy_points=data
            
    cdef double intersect_c(self, vector_t p1, vector_t p2):
        cdef:
            double max_length = sep_(p1, p2)
            double h = (self.z_plane-p1.z)/(p2.z-p1.z)
            double X, Y
        
        if (h<self.tolerance) or (h>1.0):
            #print "H", h
            return 0
        X = p1.x + h*(p2.x-p1.x)
        Y = p1.y + h*(p2.y-p1.y)
        #test for (X,Y) in polygon
        if point_in_polygon_c(X,Y, self._xy_points):
            return h * max_length
        else:
            return 0.0
        
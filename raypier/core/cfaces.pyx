#!/usr/bin/env python

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

"""
Cython module for Face definitions
"""

#maybe this is a cython .15 thing?
#from libc.math import INFINITY, M_PI, sqrt, pow, fabs, cos, sin, acos, atan2

cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)

cdef extern from "float.h":
    double DBL_MAX

cdef double INF=(DBL_MAX+DBL_MAX)

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, Shape, Distortion

import numpy as np
cimport numpy as np_
cimport cython
from cython.parallel import prange

cdef struct flatvector_t:
    double x,y
    
    
cdef class ShapedFace(Face):
    cdef:
        public Shape shape
        public int invert_normals
        
    def __cinit__(self, **kwds):
        self.shape = kwds.get("shape", Shape())
        self.invert_normals = int(kwds.get('invert_normals', 0))
        
    cdef double eval_z_c(self, double x, double y) nogil:
        return 0.0
    
    cdef double eval_implicit_c(self, double x, double y, double z) nogil:
        return z - self.eval_z_c(x,y)
    
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    def eval_z_extent(self, double[:] x, double[:] y):
        """Samples the face surface z-values over the 2d grid given by the x- and y-
           arrays. Returns the max and min values of z"""
        cdef:
            size_t i,j,ni=x.shape[0], nj=y.shape[0]
            double maxv, minv, v
        maxv = self.eval_z_c(x[0],y[0])+0.1
        minv = maxv-0.2
        with nogil:
            for i in range(ni):
                for j in range(nj):
                    v = self.eval_z_c(x[i],y[j])
                    if v < minv:
                        minv = v
                    if v > maxv:
                        maxv = v
        return (minv, maxv)
    
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    def eval_z_points(self, double[:,:] points):
        cdef:
            size_t i, ni=points.shape[0]
            np_.ndarray aout = np.empty((ni,), dtype='d')
            double[:] out = aout
            
        with nogil:
            ###Not using prange here because ZernikeDistortion uses a cache for its calculation. 
            ###Threading screws this up.
            for i in prange(ni):
                out[i] = self.eval_z_c(points[i,0],points[i,1])
        return aout
    
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    def eval_implicit_grid(self, double[:] x, double[:] y, double[:] z):
        cdef:
            size_t i,j,k, ni=x.shape[0], nj=y.shape[0], nk=z.shape[0]
            np_.ndarray aout = np.empty((ni,nj,nk), dtype='d')
            double[:,:,:] out = aout
            int sign
            
        if self.invert_normals:
            sign = -1
        else:
            sign = 1
            
        with nogil:
            for i in prange(ni):
                for j in range(nj):
                    for k in range(nk):
                        out[i,j,k] = sign*self.eval_implicit_c(x[i], y[j], z[k])
        return aout
    
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    def eval_implicit_points(self, double[:,:] points):
        cdef:
            size_t i, ni=points.shape[0]
            np_.ndarray aout = np.empty((ni,), dtype='d')
            double[:] out = aout
            int sign
        if self.invert_normals:
            sign = -1
        else:
            sign = 1
            
        with nogil:
            for i in prange(ni):
                out[i] = sign*self.eval_implicit_c(points[i,0],points[i,1],points[i,2])
        return aout
    

cdef class CircularFace(Face):
    cdef public double diameter, offset, z_plane
    
    params = ['diameter', 'offset']
    
    def __cinit__(self, **kwds):
        self.z_plane = kwds.get('z_plane', 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
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
        if is_base_ray and (X*X + Y*Y) > (d*d/4):
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
    
    
cdef class ShapedPlanarFace(ShapedFace):
    cdef:
        public double z_height
        
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
        """
        cdef:
            double max_length = sep_(p1, p2)
            double h = (self.z_height-p1.z)/(p2.z-p1.z)
            double X, Y
                    
        if (h<self.tolerance) or (h>1.0):
            return -1
        
        X = p1.x + h*(p2.x-p1.x)
        Y = p1.y + h*(p2.y-p1.y)
        if is_base_ray and not (<Shape>self.shape).point_inside_c(X,Y):
            return -1
        return h*max_length
        
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef vector_t normal
        normal.x=0
        normal.y=0
        normal.z=1
        return normal
    
    cdef double eval_z_c(self, double x, double y) nogil:
        return self.z_height
    
    cdef double eval_implicit_c(self, double x, double y, double z) nogil:
        return (z - self.z_height)
    
    
cdef class ElipticalPlaneFace(Face):
    cdef public double g_x, g_y, diameter
    
    params = ['diameter']
    
    def __cinit__(self, **kwds):
        self.g_x = kwds.get('g_x', 0.0)
        self.g_y = kwds.get('g_y', 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            double max_length = sep_(p1, p2)
            double h = (self.g_x*p1.x + self.g_y*p1.y - p1.z) / \
                        ((p2.z-p1.z) - self.g_x*(p2.x-p1.x) - self.g_y*(p2.y-p1.y))
            double X,Y,Z, d=self.diameter
            
        if (h<self.tolerance) or (h>1.0):
            return 0
        
        X = p1.x + h*(p2.x-p1.x)
        Y = p1.y + h*(p2.y-p1.y)
        Z = p1.z + h*(p2.z-p1.z)
        if is_base_ray and (X*X + Y*Y) > (d*d/4):
            #print "X", X, "Y", Y
            return 0
        return h * max_length
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef vector_t normal
        normal.x=self.g_x
        normal.y=self.g_y
        normal.z=-1
        return norm_(normal)
    
    
cdef class RectangularFace(Face):
    cdef public double length, width, offset, z_plane
    
    params = ['length', 'width', 'offset']
    
    def __cinit__(self, **kwds):
        self.z_plane = kwds.get('z_plane', 0.0)
        self.width = kwds.get("width", 2.0)
        self.length = kwds.get("length", 5.0)
        self.offset = kwds.get("offset", 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
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
            double X, Y, lngth=self.length, wdth = self.width
            
        if (h<self.tolerance) or (h>1.0):
            return 0
        
        if is_base_ray:
            X = p1.x + h*(p2.x-p1.x) - self.offset
            Y = p1.y + h*(p2.y-p1.y)
            
            #if x or y displacement is greater than length or width of rectangle, no intersect
            if X*X > lngth*lngth/4:
                return 0
            if Y*Y > wdth*wdth/4:
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
    
    #Don't want curvature in this list in case the owner defines the 
    #curvature in a different way
    params = ['diameter',] 
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
    
    cdef double intersect_c(self, vector_t r, vector_t p2, int is_base_ray):
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
            double A,B,C,D, cz, a1, a2
            vector_t s, d, pt1, pt2
            
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
        
        #1st root
        a1 = (-B+D)/(2*A) 
        pt1 = addvv_(r, multvs_(s, a1))
        #2nd root
        a2 = (-B-D)/(2*A)
        pt2 = addvv_(r, multvs_(s, a2))
        
        if self.curvature >= 0:
            if pt1.z < cz:
                a1 = INF
            if pt2.z < cz:
                a2 = INF
        else:
            if pt1.z > cz:
                a1 = INF
            if pt2.z > cz:
                a2 = INF
            
        D = self.diameter*self.diameter/4.
        
        if is_base_ray:
            if (pt1.x*pt1.x + pt1.y*pt1.y) > D:
                a1 = INF
            if (pt2.x*pt2.x + pt2.y*pt2.y) > D:
                a2 = INF
        
        if a2 < a1:
            a1 = a2
        
        if a1>1.0 or a1<self.tolerance:
            return 0
        return a1 * sep_(r, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        
        p.z -= (self.z_height - self.curvature)
        if self.curvature < 0:
            p.z = -p.z
            p.y = -p.y
            p.x = -p.x
        return norm_(p)
    
    
cdef class ShapedSphericalFace(ShapedFace):
    cdef:
        public double curvature, z_height
    
    #Don't want curvature in this list in case the owner defines the 
    #curvature in a different way
    params = [] 
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
        self.curvature = kwds.get("curvature", 100.0)
    
    cdef double intersect_c(self, vector_t r, vector_t p2, int is_base_ray):
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
            double A,B,C,D, cz, a1, a2
            vector_t s, d, pt1, pt2
            
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
        
        #1st root
        a1 = (-B+D)/(2*A) 
        pt1 = addvv_(r, multvs_(s, a1))
        #2nd root
        a2 = (-B-D)/(2*A)
        pt2 = addvv_(r, multvs_(s, a2))
        
        if self.curvature >= 0:
            if pt1.z < cz:
                a1 = INF
            if pt2.z < cz:
                a2 = INF
        else:
            if pt1.z > cz:
                a1 = INF
            if pt2.z > cz:
                a2 = INF
           
        if is_base_ray:
            if not (<Shape>self.shape).point_inside_c( pt1.x, pt1.y ):
                a1 = INF
            if not (<Shape>self.shape).point_inside_c( pt2.x, pt2.y ):
                a2 = INF
        
        if a2 < a1:
            a1 = a2
        
        if a1>1.0 or a1<self.tolerance:
            return 0
        return a1 * sep_(r, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        
        p.z -= (self.z_height - self.curvature)
        if self.curvature < 0:
            p.z = -p.z
            p.y = -p.y
            p.x = -p.x
        return norm_(p)
    
    cdef double eval_z_c(self, double x, double y) nogil:
        cdef:
            double r2 = x*x + y*y
            double c = self.curvature
        if c >= 0:
            return self.z_height + sqrt(c*c - r2) - c
        else:
            return self.z_height - sqrt(c*c - r2) - c
    
    cdef double eval_implicit_c(self, double x, double y, double z) nogil:
        cdef:
            double cz, c=self.curvature

        cz = z - self.z_height + c
        
        if c >= 0:
            return ((x*x + y*y + cz*cz) - (c*c))
        else:
            return -((x*x + y*y + cz*cz) - (c*c))
    
    
cdef class ExtrudedPlanarFace(Face):
    cdef:
        double x1_, y1_, x2_, y2_
        public double z1, z2
        vector_t normal
        
    def __cinit__(self, **kwds):
        self.x1 = kwds.get('x1',0)
        self.y1 = kwds.get('y1',0)
        self.x2 = kwds.get('x2',0)
        self.y2 = kwds.get('y2',0)
        self.z1 = kwds.get('z1',0)
        self.z2 = kwds.get('z2',0)
        
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
            
        n.y = self.x2_ - self.x1_
        n.x = self.y1_ - self.y2_
        n.z = 0
        self.normal = norm_(n)
        
    cdef double intersect_c(self, vector_t r, vector_t p2, int is_base_ray):
        cdef: 
            vector_t s, u, v
            double a, dz
            
        u.x = self.x1_
        u.y = self.y1_
        
        v.x = self.x2_ - u.x
        v.y = self.y2_ - u.y
        
        s = subvv_(p2, r)
        
        if is_base_ray:
            #fractional distance of intersection along edge
            a = (s.y*(u.x-r.x) - s.x*(u.y-r.y)) / (s.x*v.y - s.y*v.x)
            
            if a<0:
                return 0
            if a>1:
                return 0
        #distance of intersection along ray (in XY plane)
        a = (v.x*(r.y-u.y) - v.y*(r.x-u.x)) / (s.x*v.y - s.y*v.x)
        
        if is_base_ray:
            #distance in 3D
            dz = a*(p2.z - r.z)
            if self.z1 < (r.z+dz) < self.z2:
                return a * mag_(s)
            else:
                return 0
        else:
            return a * mag_(s)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        return self.normal

#
# Functions used for Bezier math.  should go in utils.pyx when there is one
#


    
cdef double eval_bezier(double t, double cp0, double cp1, double cp2, double cp3):
    #just evaluate a cubic bezier spline
    return cp0*((1-t)**3) + 3*cp1*t*((1-t)**2) + 3*cp2*(1-t)*(t**2) + cp3*(t**3)

cdef double dif_bezier(double t, double cp0, double cp1, double cp2, double cp3):
    #calc the derivative of a cubic bezier when parameter = t
    cdef long double A, B, C     #just doin this old school polynomial style
    A = cp3-3*cp2+3*cp1-cp0
    B = 3*cp2-6*cp1+3*cp0
    C = 3*cp1-3*cp0
    return 3*A*t**2 + 2*B*t + C

cdef struct poly_roots:
    #my apologies for highly unusable code, 
    #I am such a noob at passing values between functions
    double roots[3] 
    int n
    
cdef poly_roots roots_of_cubic(double a, double b, double c, double d):
    #this code is known not to work in the case of (x-c)^3 (triple zero)
    # **TODO ** fix this
    # TODO: cubic solution explodes with small a, co quadratic is used.
    # using newton's root finding method could quickly polish up answer.
    cdef:
        long double a1 = b/a, a2 = c/a, a3 = d/a
        long double Q = (a1*a1 - 3.0*a2)/9.0
        long double R = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3)/54.0
        long double R2_Q3 = R*R - Q*Q*Q
        long double theta
        poly_roots x= ((0.0,0.0,0.0),0)
    if fabs(a) <= 0.0000000001:
        #^this, precision is less than ideal here
        if fabs(b) <= 0.0000000001:
            if c == 0:
                #our "polynomial" is a constant
                # return a nonsense answer
                x.n = 1
                x.roots[0] = 0
            else:
                #print "solved for line"
                #our polynomial is a line cx+d
                x.n = 1
                x.roots[0] = -d/c
        else:
            #print "solved for quadratic"
            #our ploynomial is quadratic, bx**2+cx+d
            a1 = c**2 - 4*b*d
            a1 = sqrt(a1)
            x.n = 2
            x.roots[0] = (-c+a1)/(2*b)
            x.roots[1] = (-c-a1)/(2*b)
            #print "roots: ",x.roots[0],x.roots[1]
    else:
        #our polynomial is actually cubic   
        #original code had <= here, but if determinant is 0 then only one answer
        if (R2_Q3 < 0):
            x.n = 3
            theta = acos(R/sqrt(Q*Q*Q))
            x.roots[0] = -2.0*sqrt(Q)*cos(theta/3.0) - a1/3.0;
            x.roots[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a1/3.0
            x.roots[2] = -2.0*sqrt(Q)*cos((theta+4.0*M_PI)/3.0) - a1/3.0
            #print "triple",x.roots[0],x.roots[1],x.roots[2]
        else:
            x.n = 1
            a2 = pow(sqrt(R2_Q3)+fabs(R), 1/3.0)
            a2 += Q/x.roots[0]
            a2 *= (1 if (R < 0.0) else -1)
            a2 -= a1/3.0
            x.roots[0] = a2 
            #print "single: ",x.roots[0]
    return x

cdef flatvector_t rotate2D(double phi, flatvector_t p):
    cdef flatvector_t result
    result.x = p.x*cos(phi) - p.y*sin(phi)
    result.y = p.x*sin(phi) + p.y*cos(phi)     
    return result

cdef class ExtrudedBezierFace(Face):

    cdef:
        double z_height_1, z_height_2
        flatvector_t mincorner, maxcorner         #corners of x-y box that bounds entire spline
        np_.ndarray curves_array        

    cdef int ccw(self, flatvector_t  A, flatvector_t  B, flatvector_t  C):
        #used by an ingenious line segment intersect algorithim I found.
        #determines counter clockwiseness of points
        return (C.y-A.y)*(B.x-A.x) > (B.y-A.y)*(C.x-A.x)

    cdef int line_seg_overlap(self, flatvector_t A, flatvector_t B, flatvector_t C, flatvector_t D):
        #check if two line segments overlap eachother.  Used for rough tests of
        #intersection before committing to much computation to the potential intersection.
        # A and B are the begin and end of one line; C,D the other. 
        # Code from http://www.bryceboe.com/2006/10/23/line-segment-intersection-algorithm/
        #AB crosses CD if ABCD are all cw or ccw.
        return self.ccw(A,C,D) != self.ccw(B,C,D) and self.ccw(A,B,C) != self.ccw(A,B,D)
        
    cdef int pnt_in_hull(self,flatvector_t p, flatvector_t A, flatvector_t B, flatvector_t C, flatvector_t D):
        #instead of convex hull, just use xy bounding box, which is quicker to construct

        cdef int i,j,k
        cdef float w

        #if p is bigger than atleast one but not all, the it is in the box
        i = p.x > A.x or p.x>B.x or p.x>C.x or p.x>D.x
        j = p.x > A.x and p.x>B.x and p.x>C.x and p.x>D.x
        k = i and not j
        i = p.y > A.y or p.y>B.y or p.y>C.y or p.y>D.y
        j = p.y > A.y and p.y>B.y and p.y>C.y and p.y>D.y
        i = i and not j

        #maybe a bad idea, but this is incase bezier is actually a plane in x or y only:
        w = A.x - D.x
        w = w*w
        if w <= .0005:
            w = B.x - A.x
            w=w*w
            if w <= .0005:
                i=k=True
        else:
            w = A.y - D.y
            w = w*w
            if w <= .0005:
                w = B.y - A.y
                w=w*w
                if w <= .0005:
                    i=k=True
        return i and k

    def __cinit__(self,np_.ndarray[np_.float64_t ,ndim=3] beziercurves,double z_height_1=0,double z_height_2=0, **kwds):
        cdef flatvector_t temp1, temp2
        self.curves_array = beziercurves
        self.z_height_1 = z_height_1
        self.z_height_2 = z_height_2
        temp1.x = beziercurves[0,0,0]
        temp1.y = beziercurves[0,0,1]
        temp2.x = beziercurves[0,0,0]
        temp2.y = beziercurves[0,0,1]

        for bezierpts in beziercurves:
            for pair in bezierpts:
                if pair[0] < temp1.x: temp1.x=pair[0]
                if pair[0] > temp2.x: temp2.x=pair[0]
                if pair[1] < temp1.y: temp1.y=pair[1]
                if pair[1] > temp2.y: temp2.y=pair[1]

        self.mincorner = temp1
        self.maxcorner = temp2

    cdef double intersect_c(self, vector_t ar, vector_t pee2, int is_base_ray):

        cdef: 
            flatvector_t tempvector
            flatvector_t r, p2, s, origin
            flatvector_t cp0,cp1,cp2,cp3  #holds control points for spline segment under scrutiny
            double result = INF            #length of ray before it intersects surface. 0 if no valid intersection
            double dZ                       #rate of change of z. dZ*result+Z0 gives Z coordinate
            double A,B,C,D,t,a,b,c,d
            poly_roots ts
            vector_t    tempv
        #print "\ncalled intersection"

        origin.x = 0 
        origin.y = 0
        #first off, does ray even enter the depth of the extrusion?
        if (ar.z < self.z_height_1 and pee2.z <self.z_height_1) or (ar.z > self.z_height_2 and pee2.z > self.z_height_2):
            return 0    #does not
        #print "ray passes z test"
        #strip useless thrid dimension from ray vector
        r.x = ar.x
        r.y = ar.y
        p2.x = pee2.x
        p2.y = pee2.y


        #check if ray intersects x-y bounding box of spline at all
        tempvector.x = self.mincorner.x
        tempvector.y = self.maxcorner.y
        if not self.line_seg_overlap(r,p2,self.mincorner,tempvector):
            if not  self.line_seg_overlap(r,p2,tempvector,self.maxcorner):
                tempvector.x = self.maxcorner.x
                tempvector.y = self.mincorner.y
                if not self.line_seg_overlap(r,p2,self.maxcorner,tempvector):
                    if not self.line_seg_overlap(r,p2,tempvector,self.mincorner):
                        return 0    #no intersections
        #print "ray is in big bounding box"
        #segment intersects with gross bounding box,
        #Calc dZ and the 2D origin adjusted ray, because they will probably be used.
        tempv = subvv_ (pee2,ar) 
        dZ = tempv.z
        s.x=tempv.x
        s.y=tempv.y
        theta = atan2(s.y,s.x)
                
        s = rotate2D(-theta,s)
        # now, loop through curves and see if segment 1) intersects with individual convex hulls
        # 2) intersects with spline (return points)
        for curve in self.curves_array.copy():            #load up control points
            for pt in curve:
                pt[0]=pt[0]-ar.x
                pt[1]=pt[1]-ar.y 
            cp0.x,cp0.y = curve[0].copy()
            cp1.x,cp1.y = curve[1].copy()
            cp2.x,cp2.y = curve[2].copy()
            cp3.x,cp3.y = curve[3].copy()
            
            
            #rotate ctrl points such that ray is along the x axis
            cp0 = rotate2D(-theta,cp0)
            cp1 = rotate2D(-theta,cp1)
            cp2 = rotate2D(-theta,cp2)
            cp3 = rotate2D(-theta,cp3)
            #test for intersection between ray (actually segment) and convex hull
            if self.line_seg_overlap(origin,s,cp0,cp1) or self.line_seg_overlap(origin,s,cp1,cp2) or self.line_seg_overlap(origin,s,cp2,cp3) or self.line_seg_overlap(origin,s,cp3,cp0):
                #print "inside intersect hull"
                #Ray does intersect this convex hull.  Find solution:
                #Setup A,B,C and D (bernstein polynomials)
                A = cp3.y-3*cp2.y+3*cp1.y-cp0.y
                B = 3*cp2.y-6*cp1.y+3*cp0.y
                C = 3*cp1.y-3*cp0.y
                D = cp0.y
                #solve for t
                #print "ABCD: ",A,B,C,D
                ts = roots_of_cubic(A,B,C,D)
                while ts.n > 0:
                    ts.n-=1
                    t = ts.roots[ts.n]
                    #print "root: ", t
                    #make sure solution is on valid interval
                    if 0.<t<1.:
                        #the x value will also be the length, which is the form of result
                        b = eval_bezier(t,cp0.x,cp1.x,cp2.x,cp3.x)
                        #print "b at t",b,t
                        #is x within bounds?
                        if 0 < b < s.x:
                            #print "in range"
                            #is point within Z bounds?
                            c = dZ*b/s.x   
                            a = c+ar.z     
                            if self.z_height_1 < a < self.z_height_2:
                                #print "in z: ",B,result
                                #is this the shortest length to an intersection so far?
                                b = sqrt(c**2+b**2)
                                #print "this: ",b
                                if b < result and b > self.tolerance:
                                    result = b
                                    #print "b at t",b,t

        if result == INF: result = 0

        return result



    cdef vector_t compute_normal_c(self, vector_t p):
        cdef:
            flatvector_t ray,cp0,cp1,cp2,cp3,rotated
            double theta, tmp, t
            poly_roots ts
            np_.ndarray box
        #print "called looking for normal",p.x,p.y
        ray.x = p.x
        ray.y = p.y
        theta = atan2(p.y,p.x)
        #find which curve this point is in
        for curve in self.curves_array.copy():            #load up control points
            cp0.x,cp0.y = curve[0].copy()
            cp1.x,cp1.y = curve[1].copy()
            cp2.x,cp2.y = curve[2].copy()
            cp3.x,cp3.y = curve[3].copy()

            #is point even in this hull?
            #print "tried at least"
            if self.pnt_in_hull(ray,cp0,cp1,cp2,cp3):
                #then, solve for t
                #print "in hull"
                cp0 = rotate2D(-theta,cp0)
                cp1 = rotate2D(-theta,cp1)
                cp2 = rotate2D(-theta,cp2)
                cp3 = rotate2D(-theta,cp3)
                
                #Setup A,B,C and D (bernstein polynomials)
                A = cp3.y-3*cp2.y+3*cp1.y-cp0.y
                B = 3*cp2.y-6*cp1.y+3*cp0.y
                C = 3*cp1.y-3*cp0.y
                D = cp0.y
                ts = roots_of_cubic(A,B,C,D)
                
                #print "normal roots: ",ts.n,ts.roots[0],ts.roots[1],ts.roots[2]
                while ts.n > 0:
                    ts.n -=1
                    t = ts.roots[ts.n]
                    #make sure solution is within interval
                    #print "normal t: ",t
                    if 0<=t<=1:
                        #ok, then is this the t to the same point p? 
                        tmp = eval_bezier(t,cp0.x,cp1.x,cp2.x,cp3.x)
                        #print "normal b at t: ",tmp,t
                        #I will generously allow for rounding error 
                        #print "got here with point: ",tmp,
                        if tmp**2 - (ray.x**2+ray.y**2) < .0001:
                            #this is the single solution. return the derivative dy/dx = dy/dt / dx/dt
                            #print "that was it!"
                            ray.x = dif_bezier(t,cp0.x,cp1.x,cp2.x,cp3.x)
                            ray.y = dif_bezier(t,cp0.y,cp1.y,cp2.y,cp3.y)
                            ray = rotate2D(theta,ray)
                            p.z = 0     #trough has no slope in z
                            #direction of normal is to the left of the parametric curve
                            #slope of normal is -dx/dy
                            
                            if ray.y==0:
                                p.x=0
                                p.y = (1 if ray.x>0 else -1)
                            elif ray.y>0:
                                p.x=-1
                                p.y = ray.x/ray.y
                            elif ray.y<0:
                                p.x = 1
                                p.y = -ray.x/ray.y    
                            return norm_(p)


        #how did you get here?  p was supposed to be a point on the curve!
        print("error: Bezier normal not found, point not actually on curve!")
        p.x=p.y=p.z = 0
        return p


    
cdef int point_in_polygon_c(double X, double Y, object obj):
    cdef int i, size, ct=0
    cdef double y1, y2, h, x, x1, x2
    cdef np_.ndarray[np_.float64_t, ndim=2] pts=obj
    
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
                ct = not ct
        y1 = y2
        x1 = x2
    return ct


def point_in_polygon(double X, double Y, point_list):
    pts = np.ascontiguousarray(point_list, dtype=np.float64)
    assert pts.shape[1]==2
    return bool(point_in_polygon_c(X, Y, pts))
    
    
cdef class PolygonFace(Face):
    cdef public double z_plane
    cdef object _xy_points
    
    def __cinit__(self, z_plane=0.0, xy_points=[[]], **kwds):
        self.z_plane = z_plane
        self.xy_points = xy_points
    
    property xy_points:
        def __get__(self):
            return self._xy_points
        
        def __set__(self, pts):
            data = np.ascontiguousarray(pts, dtype=np.float64).reshape(-1,2)
            self._xy_points=data
            
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
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
        if is_base_ray and point_in_polygon_c(X,Y, self._xy_points)==1:
            return h * max_length
        else:
            return 0.0
        
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef vector_t normal
        normal.x=0
        normal.y=0
        normal.z=-1
        return normal
    
    
cdef class OffAxisParabolicFace(Face):
    cdef:
        public double EFL, diameter, height
                
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
        """
        cdef:
            double max_length = sep_(p1, p2)
            double A = 1 / (2*self.EFL), efl = self.EFL
            double a,b,c,d
            vector_t s, r, pt1, pt2
            
            
        s = subvv_(p2,p1)
        r = p1
        r.z += efl/2.
                
        a = A*(s.x**2 + s.y**2)
        b = 2*A*(r.x*s.x + r.y*s.y) - s.z
        c = A*(r.x**2 + r.y**2) - r.z
        d = b**2 - 4*a*c
        
        ###FIXME
        if d<0: #no intersection ###BROKEN
            return 0
        
        if a < 1e-10: #approximate to zero if we're close to the parabolic axis
            a1 = -c/b
            pt1 = addvv_(r, multvs_(s, a1))
            pt1.x -= efl
            if (pt1.x*pt1.x + pt1.y*pt1.y) > (self.diameter/2):
                return 0
            if a1>1.0 or a1<self.tolerance:
                return 0
            return a1 * sep_(p1, p2)
        
        else:
            d = sqrt(d)
           
            #1st root
            a1 = (-b+d)/(2*a) 
            pt1 = addvv_(r, multvs_(s, a1))
            #2nd root
            a2 = (-b-d)/(2*a)
            pt2 = addvv_(r, multvs_(s, a2))
                           
            pt1.x -= efl
            pt2.x -= efl
            
            if is_base_ray:
                d = self.diameter
                d *= d/4.
                
                if (pt1.x*pt1.x + pt1.y*pt1.y) > d:
                    a1 = INF
                if (pt2.x*pt2.x + pt2.y*pt2.y) > d:
                    a2 = INF
            
            if a2 < a1:
                a1 = a2
            
            if a1>1.0 or a1<self.tolerance:
                return 0
            return a1 * sep_(p1, p2)
#
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            vector_t normal
            double A = 1 / (2*self.EFL)
            double B, dz, m2
        
        m2 = p.x*p.x + p.y*p.y
        B = 4*m2*A*A
        dz = -sqrt( B/(B+1) )
        m2 = sqrt(m2)
        
        normal.x= -(dz * p.x)/m2
        normal.y= -(dz * p.y)/m2
        normal.z= -1 / sqrt( B+1 )
        return normal
    
        
cdef class EllipsoidalFace(Face):
    cdef:
        public double major, minor #axis lengths
        transform_t trans, inv_trans
        public double x1, x2, y1, y2, z1, z2 #local bounds of the ellpsoid block
        
    property transform:
        def __get__(self):
            t = Transform()
            t.trans = self.trans
            return t
        
        def __set__(self, Transform t):
            self.trans = t.trans
            
    property inverse_transform:
        def __set__(self, Transform t):
            self.inv_trans = t.trans
           
        def __get__(self):
            cdef Transform t=Transform()
            t.trans = self.inv_trans
            return t
            
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            double B,A, a, b, c, d
            
            vector_t S = subvv_(p2, p1)
            vector_t r = transform_c(self.trans, p1)
            vector_t s = transform_c(self.trans, p2)
            
        s = subvv_(s, r)
        
        B = self.minor**2
        A = self.major**2
        
        a = A*(s.z*s.z + s.y*s.y) + B*s.x*s.x
        b = 2*( A*(r.z*s.z + r.y*s.y) + B*r.x*s.x )
        c = A*(r.z*r.z + r.y*r.y) + B*r.x*r.x - A*B
        
        d = b*b - 4*a*c
        d = sqrt(d)
        root1 = (-b + d)/(2*a)
        root2 = (-b - d)/(2*a)
        p2 = addvv_(p1, multvs_(S, root2))
        p1 = addvv_(p1, multvs_(S, root1))
        
        if is_base_ray:
            if not self.x1 < p2.x < self.x2:
                root2 = 2
            if not self.y1 < p2.y < self.y2:
                root2 = 2
            if not self.z1 < p2.z < self.z2:
                root2 = 2
                
            if not self.x1 < p1.x < self.x2:
                root1 = 2
            if not self.y1 < p1.y < self.y2:
                root1 = 2
            if not self.z1 < p1.z < self.z2:
                root1 = 2
        
        if root1 < self.tolerance:
            root1 = 2
        if root2 < self.tolerance:
            root2 = 2
        if root1 > root2:
            root1 = root2
        if root1 > 1:
            return 0
        return root1*mag_(S)
        
    cdef vector_t compute_normal_c(self, vector_t p):
        cdef vector_t n
        
        p = transform_c(self.trans, p)
        
        n.x = p.x/-(self.major**2)
        n.y = p.y/-(self.minor**2)
        n.z = p.z/-(self.minor**2)
        
        n = rotate_c(self.inv_trans, n)
        return norm_(n)
    
    def update(self):
        super(EllipsoidalFace, self).update()
        owner = self.owner
        self.sync_transform(owner.ellipse_trans)
        self.major, self.minor = owner.axes
        self.x1, self.x2 = owner.X_bounds
        self.y1, self.y2 = owner.Y_bounds
        self.z1, self.z2 = owner.Z_bounds
        
    def sync_transform(self, vtk_trans):
        m = vtk_trans.matrix
        rot = [[m.get_element(i,j) for j in xrange(3)] for i in xrange(3)]
        dt = [m.get_element(i,3) for i in xrange(3)]
        #print "TRANS", rot, dt
        self.transform = Transform(rotation=rot, translation=dt)
        inv_trans = vtk_trans.linear_inverse
        m = inv_trans.matrix
        rot = [[m.get_element(i,j) for j in xrange(3)] for i in xrange(3)]
        dt = [m.get_element(i,3) for i in xrange(3)]
        self.inverse_transform = Transform(rotation=rot, translation=dt)
        
        
cdef class SaddleFace(ShapedFace):
    cdef:
        public double z_height, curvature
        
    params=[]
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get("z_height", 0.0)
        self.curvature = kwds.get("curvature", 0.0)
        
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            double A=sqrt(6.0), root, denom, a1, a2
            vector_t d,p, pt1, pt2
            
        A *= self.curvature
        p = p1
        p.z -= self.z_height
        d = subvv_(p2,p1)
        
        if d.x==0.0:
            a1 = (-A*(p.x*p.y) + p.z)/(A*d.y*p.x - d.z)
            a2 = INF
        elif d.y==0.0:
            a1 = (-A*(p.x*p.y) + p.z)/(A*d.x*p.y - d.z)
            a2 = INF
        else:
            root = (A**2)*(d.x**2)*(p.y**2) - 2*(A**2)*d.x*d.y*p.x*p.y + (A**2)*(d.y**2)*(p.x**2) + \
                    4*A*d.x*d.y*p.z - 2*A*d.x*d.z*p.y - 2*A*d.y*d.z*p.x + d.z**2
                    
            if root < 0:
                return -1
            
            root = sqrt(root)
            denom = 2*A*(d.x*d.y)
            
            a1 = a2 = -A*d.x*p.y - A*d.y*p.x + d.z
            a1 += root
            a2 -= root
            a1 /= denom
            a2 /= denom
        
        #print("a1:", a1, "a2:", a2, "root:", root, "sep:", sep_(p1,p2))
        
        pt1 = addvv_(p1, multvs_(d, a1))
        pt2 = addvv_(p1, multvs_(d, a2))
        
        if a1 < 0.0:
            a1 = INF
        if a2 < 0.0:
            a2 = INF
        
        if is_base_ray:
            if not (<Shape>self.shape).point_inside_c( pt1.x, pt1.y ):
                a1 = INF
            if not (<Shape>self.shape).point_inside_c( pt2.x, pt2.y ):
                a2 = INF
        
        if a2 < a1:
            a1 = a2
        
        if a1>1.0 or a1<self.tolerance:
            return -1
        return a1 * sep_(p1, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            vector_t n
            double rt6 = sqrt(6)*self.curvature
        
        n.z = 1.0
        n.x = -rt6*p.y
        n.y = -rt6*p.x
        return norm_(n)
    
    cdef double eval_z_c(self, double x, double y) nogil:
        return sqrt(6) * self.curvature * x * y
        
        
    
        
cdef class CylindericalFace(ShapedFace):
    cdef:
        public double z_height, radius
    
    params = [] 
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
        self.radius = kwds.get("radius", 100.0)
        
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            double a1, a2, cz, ox2, oz2, dx2, dz2, denom, R=self.radius
            double R2 = R*R
            vector_t d,o, pt1, pt2
            
        o = p1
        o.z -= self.z_height
        d = subvv_(p2,p1)
        
        ox2 = o.x*o.x
        oz2 = o.z*o.z
        dx2 = d.x*d.x
        dz2 = d.z*d.z
        
        root = R2*dz2 - 2*R*dx2*o.z + 2*R*d.x*d.z*o.x - dx2*oz2 + 2*d.x*d.z*o.x*o.z - dz2*ox2
        
        if root < 0:
            return -1
        
        root = sqrt(root)
        denom = dx2 + dz2
        
        a1 = a2 = -R*d.z - d.x*o.x - d.z*o.z
        a1 += root
        a2 -= root
        a1 /= denom
        a2 /= denom
        
        pt1 = addvv_(p1, multvs_(d, a1))
        pt2 = addvv_(p1, multvs_(d, a2))
        
        cz = self.z_height - self.radius
        
        if R >= 0:
            if pt1.z < cz:
                a1 = INF
            if pt2.z < cz:
                a2 = INF
        else:
            if pt1.z > cz:
                a1 = INF
            if pt2.z > cz:
                a2 = INF
        
        if is_base_ray:
            if not (<Shape>self.shape).point_inside_c( pt1.x, pt1.y ):
                a1 = INF
            if not (<Shape>self.shape).point_inside_c( pt2.x, pt2.y ):
                a2 = INF
        
        if a2 < a1:
            a1 = a2
        
        if a1>1.0 or a1<self.tolerance:
            return -1
        return a1 * sep_(p1, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        
        p.z -= (self.z_height - self.radius)
        if self.radius < 0:
            p.z = -p.z
            p.x = -p.x
        p.y = 0
        return norm_(p)
    
    cdef double eval_z_c(self, double x, double y) nogil:
        cdef:
            double R = self.radius
        if R>=0:
            return self.z_height + sqrt(R*R - x*x) - R
        else:
            return self.z_height - sqrt(R*R - x*x) - R
    
    
cdef class AxiconFace(ShapedFace):
    """
    While technically, we can use the conic surface to generate a cone, it requires setting some parameters to infinity which 
    is often inaccurate to compute.
    
    The gradient is the slope of the sides, dz/dr
    """
    cdef:
        public double z_height, gradient
        
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
        self.gradient = kwds.get('gradient', 0.0)
        
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            double a1, a2, root, ox2, oy2, oz2, dx2, dy2, dz2, beta2, denom
            double beta = self.gradient
            vector_t o, d, pt1, pt2
            
        d = subvv_(p2, p1)
        o = p1
        o.z -= self.z_height
        
        beta2 = beta*beta
        ox2 = o.x*o.x
        oy2 = o.y*o.y
        oz2 = o.z*o.z
        dx2 = d.x*d.x
        dy2 = d.y*d.y
        dz2 = d.z*d.z
        
        root = -beta2*dx2*oy2 + 2*beta2*d.x*d.y*o.x*o.y - beta2*dy2*ox2 + dx2*oz2 - 2*d.x*d.z*o.x*o.z + \
                dy2*oz2 - 2*d.y*d.z*o.y*o.z + dz2*ox2 + dz2*oy2
        denom = (beta2*dx2 + beta2*dy2 - dz2)
                
        if root < 0: #no intersection
            return -1
        
        root = beta*sqrt(root)
        
        a1 = -beta2*d.x*o.x - beta2*d.y*o.y + d.z*o.z
        a2 = a1 + root
        a1 -= root
        a1 /= denom
        a2 /= denom
        
        pt1 = addvv_(p1, multvs_(d, a1))
        pt2 = addvv_(p1, multvs_(d, a2))
        
        if pt1.z > self.z_height:
            a1 = INF
        if pt2.z > self.z_height:
            a2 = INF
        
        if is_base_ray:
            if not (<Shape>self.shape).point_inside_c( pt1.x, pt1.y ):
                a1 = INF
            if not (<Shape>self.shape).point_inside_c( pt2.x, pt2.y ):
                a2 = INF
        
        if a2 < a1:
            a1 = a2
        
        if a1>1.0 or a1<self.tolerance:
            return -1
        return a1 * sep_(p1, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            double r, beta
            
        beta = self.gradient
        r = sqrt(p.x*p.x + p.y*p.y)
        p.z = 1.0
        p.x = beta * p.x/r
        p.y = beta * p.y/r
        return p
    
    cdef double eval_z_c(self, double x, double y) nogil:
        return self.z_height - (self.gradient * sqrt(x*x + y*y))
    
    
cdef double intersect_conic(vector_t a, vector_t d, double curvature, double conic_const):
    cdef:
        double beta = 1 + conic_const
        double R = -curvature
        double A,B,C,D
        #double pt1, pt2
    
    #print("inputs:", a.x, a.y, a.z, d.x, d.y, d.z, beta, R)
    ###Obtained by sympy evaluation of the equations
    ### Arranged into quadratic form i.e. A*(alpha**2) + B*alpha + C = 0
    A = beta**2*d.z**2 + beta*d.x**2 + beta*d.y**2
    B = -2*R*beta*d.z + 2*a.x*beta*d.x + 2*a.y*beta*d.y + 2*a.z*beta**2*d.z
    C = -2*R*a.z*beta + a.x**2*beta + a.y**2*beta + a.z**2*beta**2
    
    ##Get roots
    D = B*B - 4*A*C
    if D < 0: #no solution
        return -1
    
    D = sqrt(D)
    
### This turns out to not be necessary
#     #1st root
#     a1 = (-B+D)/(2*A) 
#     pt1 = a.z + d.z*a1 #addvv_(a, multvs_(d, a1))
#     #2nd root
#     a2 = (-B-D)/(2*A)
#     pt2 = a.z + d.z*a2 #addvv_(a, multvs_(d, a2))
#     
#     print("R.beta.dz:", R*beta*d.z, a1, pt1, a2, pt2, A, B, C, D)
#     if R*beta*d.z <= 0:
#         if pt1 < R/beta:
#             a1 = INF
#         if pt2 < R/beta:
#             a2 = INF
#     else:
#         if pt1 > R/beta:
#             a1 = INF
#         if pt2 > R/beta:
#             a2 = INF
#             
#     if a2 < a1:
#         return a2
#     else:
#         return a1
    
    #print("R.beta.dz:", R*beta*d.z)
    if R*beta*d.z <= 0:
        #1st root
        return (-B+D)/(2*A) 
    else:
        #2nd root
        return (-B-D)/(2*A)
    
         
    
cdef class ConicRevolutionFace(ShapedFace):
    """This is surface of revolution formed from a conic section. Spherical and ellipsoidal faces
    are a special case of this.
    
    curvature = radius of curvature
    """
    cdef:
        public double curvature, z_height, conic_const
    
    params = [] 
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
        self.conic_const = kwds.get('conic_const', 0.0)
        self.curvature = kwds.get('curvature', 10.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
        """
        cdef:
            double a1
            vector_t d, a, pt1
            
        d = subvv_(p2, p1) #the input ray direction, in local coords.
        a = p1
        a.z -= self.z_height
        
        a1 = intersect_conic(a, d, self.curvature, self.conic_const)
        
        pt1 = addvv_(a, multvs_(d, a1))
            
        if is_base_ray and not (<Shape>(self.shape)).point_inside_c(pt1.x, pt1.y):
            return INF
            
        if a1>1.0 or a1<self.tolerance:
            return -1
        
        return a1 * sep_(p1, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            double R = -self.curvature
            double beta = 1 + self.conic_const
            vector_t g #output gradient vector
            int sign = -1 if self.invert_normals else 1
            
        p.z -= self.z_height
            
        g.z = 2*beta*(R-beta*p.z)
        g.x = - p.x * 2 * beta 
        g.y = - p.y * 2 * beta
        
        if (R*beta) < 0:
            sign *= -1
            
        g.z *= sign
        g.y *= sign
        g.x *= sign
        
        return norm_(g)
    
    cdef double eval_z_c(self, double x, double y) nogil:
        cdef:
            double r2 = (x*x) + (y*y)
            double R = self.curvature
            double R2 = R*R
        
        return self.z_height  - (r2/(R + sqrt(R2 - (1+self.conic_const)*r2)))
    
    cdef double eval_implicit_c(self, double x, double y, double z) nogil:
        return z - self.eval_z_c(x,y)
    
    
cdef struct aspheric_t:
    double R #curvature
    double beta #=1+conic const
    double A4
    double A6
    double A8
    double A10
    double A12
    double A14
    double A16
    vector_t a
    vector_t d
    
    
cdef double eval_aspheric_impf(aspheric_t A, double alpha):
    cdef: 
        double out
        double r2 = ((A.a.x + alpha*A.d.x)**2 + (A.a.y + alpha*A.d.y)**2)
    
    out = r2
    out /= A.R*(1 + sqrt(1 - A.beta*r2/(A.R**2)) )
    out -= A.a.z + alpha*A.d.z
    out += A.A4*(r2**2) + A.A6*(r2**3) + A.A8*(r2**4) + A.A10*(r2**5) + A.A12*(r2**6) + A.A14*(r2**7) + A.A16*(r2**8)
    return out


cdef double eval_aspheric_grad(aspheric_t A, double alpha):
    cdef:
        double out
        double r2 = ((A.a.x + alpha*A.d.x)**2 + (A.a.y + alpha*A.d.y)**2)
        double dx = A.d.x*(A.a.x + alpha*A.d.x)
        double dy = A.d.y*(A.a.y + alpha*A.d.y)
        
    out = A.A10*(10*dx + 10*dy) *(r2**4)
    out += A.A12*(12*dx + 12*dy) *(r2**5)
    out += A.A14*(14*dx + 14*dy) *(r2**6)
    out += A.A16*(16*dx + 16*dy) *(r2**7)
    out +=  A.A4*(4*dx + 4*dy)*(r2)  
    out +=  A.A6*(6*dx + 6*dy)*(r2**2)  
    out +=  A.A8*(8*dx + 8*dy)*(r2**3) - A.d.z  
    out += (2*dx + 2*dy)/(A.R*(sqrt(1 - A.beta*(r2)/A.R**2) + 1))  
    out += A.beta*(2*dx + 2*dy)*(r2)/(2*(A.R**3)*sqrt(1 - A.beta*r2/(A.R**2))*(sqrt(1 - A.beta*r2/(A.R**2)) + 1)**2)
    return out

    
cdef class AsphericFace(ShapedFace):
    """This is the general aspheric lens surface formula.
    
    curvature = radius of curvature
    """
    cdef:
        public double curvature, z_height, conic_const, A4, A6, A8, A10, A12, A14, A16
        public double atol
    
    params = [] 
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
        self.conic_const = kwds.get('conic_const', 0.0)
        self.curvature = kwds.get('curvature', 25.0)
        self.A4 = kwds.get('A4',0.0)
        self.A6 = kwds.get('A6',0.0)
        self.A8 = kwds.get('A8',0.0)
        self.A10 = kwds.get('A10',0.0)
        self.A12 = kwds.get('A12',0.0)
        self.A14 = kwds.get('A14',0.0)
        self.A16 = kwds.get('A16',0.0)
        self.atol = kwds.get("atol", 1.0e-8)
        
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
        """
        cdef:
            double a1, dz, f, f_last
            double tol = self.atol**2
            vector_t d, a, pt1
            aspheric_t A
            
        d = subvv_(p2, p1) #the input ray direction, in local coords.
        a = p1
        a.z -= self.z_height
        
        a1 = intersect_conic(a, d, self.curvature, self.conic_const)
        
        A.R = -self.curvature
        A.beta = 1 + self.conic_const
        A.A4 = self.A4
        A.A6 = self.A6
        A.A8 = self.A8
        A.A10 = self.A10
        A.A12 = self.A12
        A.A14 = self.A14
        A.A16 = self.A16
        A.a = a
        A.d = d

        ### Find root using Newton's method. Typically only a 3-4 iterations are required.
         
        f = f_last = eval_aspheric_impf(A, a1)
        dz = - f / eval_aspheric_grad(A, a1)
        #print("Start:", dz, a1)

        #### If we've not converged after 100 iterations, then it ain't going converge ever.
        for i in range(100):
            a1 += dz
            if dz*dz < tol:
                break
            f = eval_aspheric_impf(A, a1)
            if fabs(f) > fabs(f_last): #We're not converging
                return -1
            f_last = f
            dz = - f / eval_aspheric_grad(A, a1)
        else:
            return -1     
        
        #print("Converged:", dz, a1, i)
        
        pt1 = addvv_(a, multvs_(d, a1))

        if is_base_ray and not (<Shape>(self.shape)).point_inside_c(pt1.x, pt1.y):
            return INF
            
        if a1>1.0 or a1<self.tolerance:
            return -1
        
        #print("Ret:", a1 * sep_(p1, p2))
        return a1 * sep_(p1, p2)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            double R = -self.curvature
            double beta = 1 + self.conic_const
            vector_t g #output gradient vector
            int sign = -1 if self.invert_normals else 1
            double r2, r, df, root
            
        p.z -= self.z_height
        
        r2 = p.x*p.x + p.y*p.y
        r = sqrt(r2)
        root = sqrt(1-(beta*(r2)/(R*R)))
        df = 10*self.A10*(r2**4) + 8*self.A8*(r2**3) + 6*self.A6*(r2**2) + 4*self.A4*r2
        df += 16*self.A16*(r2**7) + 14*self.A14*(r2**6) + 12*self.A12*(r2**5) 
        df += 2/(R*(1+root))
        df += beta*(r2)/((R**3)*root*((1+root)**2))
            
        g.z = 1.0
        g.x = - df*p.x
        g.y = - df*p.y
        
        g.z *= sign
        g.y *= sign
        g.x *= sign
        
        return norm_(g)

    cdef double eval_z_c(self, double x, double y) nogil:
        cdef: 
            double out
            double r2 = (x*x) + (y*y)
            double R = self.curvature
    
        out = r2
        out /= (R + sqrt(R*R - (1+self.conic_const)*r2) )
        out += self.A4*(r2**2) + self.A6*(r2**3) + self.A8*(r2**4) + self.A10*(r2**5) + self.A12*(r2**6) + self.A14*(r2**7) + self.A16*(r2**8)
        out += self.z_height
        return out
    
    cdef double eval_implicit_c(self, double x, double y, double z) nogil:
        return z - self.eval_z_c(x,y)
    
    
cdef class DistortionFace(ShapedFace):
    """This class wraps another ShapedFace object, and applies a small distortion to it's
    surface geometry. The distortion is given by an instance of a Distortion subclass
    """
    cdef:
        public ShapedFace base_face
        public Distortion distortion
        public double accuracy
        
    def __cinit__(self, **kwds):
        face = kwds.get("base_face")
        self.base_face = face
        self.distortion = kwds.get("distortion")
        self.shape = kwds.get('shape', face.shape)
        self.accuracy = kwds.get("accuracy", 1e-6)
        
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        cdef:
            ShapedFace face=self.base_face
            double z_shift=0.0, tolerance, a1, a2, h=sep_(p2,p1)
            vector_t dxdyz, pt1, d, n, o, q1, q2
            int i
            Distortion dist=self.distortion
            Shape shape = self.shape
            
        tolerance = self.accuracy
        
        a2 = face.intersect_c(p1,p2,0) #Don't check the shape yet
        if a2>h or a2<self.tolerance: #If no intersection, then, no intersection
            #print("NO intersection", a2, p1.z, p2.z)
            return -1.0
        
        #print("start estimate:", a2)
        #starting estimate of intersection
        d = subvv_(p2,p1)
        pt1 = addvv_(p1, multvs_(d,a2/h))
        dxdyz = dist.z_offset_and_gradient_c(pt1.x, pt1.y)
        
        n = face.compute_normal_c(pt1)
        #print("base normal:", n.x, n.y, n.z)
        #print("pt1:", pt1.x, pt1.y, pt1.z)
        #print("dx:", dxdyz.x, "dy:", dxdyz.y, "z:", dxdyz.z)
        
        pt1.z += dxdyz.z
        
        n.x /= n.z
        n.y /= n.z
        n.x -= dxdyz.x
        n.y -= dxdyz.y
        
        o=subvv_(p1, pt1)
        ### Better estimate for a
        a1 = -h*dotprod_(o,n)/dotprod_(d,n)
        #print("better estimate", a1)
        if a1 < fabs(dxdyz.z):
            #print("self-intersection. bailing.")
            return -1.0
        
        for i in range(20):
            pt1 = addvv_(p1, multvs_(d,a1/h))
            #print("a_error:", a1-a2, i)
            if fabs(a1 - a2) < tolerance:
                break
            z_shift = dist.z_offset_c(pt1.x, pt1.y)
            
            q1 = p1
            q2 = p2
            q1.z -= z_shift
            q2.z -= z_shift
            a2 = face.intersect_c(q1, q2, 0)
            
            pt1 = addvv_(q1, multvs_(d,a2/h)) #this point on base
            n = face.compute_normal_c(pt1)
            dxdyz = dist.z_offset_and_gradient_c(pt1.x, pt1.y)
            pt1.z += dxdyz.z #move to distorted face
            
            n.x /= n.z
            n.y /= n.z
            n.x -= dxdyz.x
            n.y -= dxdyz.y
            
            o=subvv_(p1, pt1)
            ### Better estimate for a
            a2 = a1
            a1 = -h*dotprod_(o,n)/dotprod_(d,n)
            
        if not shape.point_inside_c(pt1.x, pt1.y):
            return -1.0
        #print("Final a1:", a1)
        return a1
            
    cdef vector_t compute_normal_c(self, vector_t p):
        cdef:
            vector_t n, p1, dxdyz

        dxdyz = self.distortion.z_offset_and_gradient_c(p.x, p.y)
        p1 = p
        p1.z -= dxdyz.z
        n = self.base_face.compute_normal_c(p1)
        n.x /= n.z
        n.y /= n.z
        n.z = 1.0
        n.x -= dxdyz.x
        n.y -= dxdyz.y
        return norm_(n)
            
    cdef double eval_z_c(self, double x, double y) nogil:
        cdef:
            double z
        z = self.base_face.eval_z_c(x,y)
        z += self.distortion.z_offset_c(x,y)
        return z
    
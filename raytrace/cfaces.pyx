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
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c

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
    
cdef class RectangularFace(Face):
    cdef public double length, width, offset, z_plane
    
    params = ['length', 'width', 'offset']
    
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
            double X, Y, lngth=self.length, wdth = self.width
            
        #print "CFACE", p1, p2
        
        if (h<self.tolerance) or (h>1.0):
            #print "H", h
            return 0
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
    
    params = ['diameter', 'curvature']
    
    def __cinit__(self, **kwds):
        self.z_height = kwds.get('z_height', 0.0)
    
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
        
        if pt1.z < cz:
            a1 = INFINITY
        if pt2.z < cz:
            a2 = INFINITY
            
        D = self.diameter*self.diameter/4.
        
        if (pt1.x*pt1.x + pt1.y*pt1.y) > D:
            a1 = INFINITY
        if (pt2.x*pt2.x + pt2.y*pt2.y) > D:
            a2 = INFINITY
        
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
        return norm_(p)
    
    
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
        #print "dist along edge:", a, r, s, v, u
        if a<0:
            return 0
        if a>1:
            return 0
        #distance of intersection along ray (in XY plane)
        a = (v.x*(r.y-u.y) - v.y*(r.x-u.x)) / (s.x*v.y - s.y*v.x)
        #print "dist along ray:", a 
        
        #distance in 3D
        dz = a*(p2.z - r.z)
        
        if self.z1 < (r.z+dz) < self.z2:
            return a * mag_(s)
        else:
            return 0
    
    cdef vector_t compute_normal_c(self, vector_t p):
        return self.normal
    
    
cdef class ExtrudedBezier2Face(Face):
    cdef:
        double X0, X1, X2, Y0, Y1, Y2, z_height_1, z_height_2
    
    cdef inline double det2x2_(self, double a11, double a12, double a21, double a22):
        cdef double out
        out = a11*a22 - a12*a21 
        return out
    
    def __cinit__(self, X0,X1,Y0,Y1,X2,Y2,z_height_1,z_height_2, **kwds):
        self.X0=X0
        self.X1=X1
        self.X2=X2
        self.Y0=Y0
        self.Y1=Y1
        self.Y2=Y2
        self.z_height_1 = z_height_1
        self.z_height_2 = z_height_2
    
    cdef double intersect_c(self, vector_t r, vector_t p2):
        cdef: 
            vector_t s
            double A, B, C, D
            double t1, t2, a1, a2
            double p0s, p1s, p2s
            double dZ, z1       #used to calculate z bounds logic
        
        dZ = p2.z-r.z
        z1 = r.z

        s = subvv_(p2, r)

        ### An inefficient implementation to begin with,
        ### to check the maths
        A = s.y*(self.X0 - 2*self.X1 + self.X2) \
            - s.x*(self.Y0 - 2*self.Y1 + self.Y2)
        B = 2*( s.x*(self.Y0-self.Y1) - s.y*(self.X0-self.X1) )
        C = s.y*(self.X0 - r.x) - s.x*(self.Y0 - r.y)
        
        D = B*B - 4*A*C
        
        if D<0: #no intersection at all
            return 0
        
        D = sqrt(D)
        
        ###two possible roots
        t1 = (-B + D) / (2*A)
        t2 = (-B - D) / (2*A)
        
        p0s = self.X0*s.x + self.Y0*s.y
        p1s = self.X1*s.x + self.Y1*s.y
        p2s = self.X2*s.x + self.Y2*s.y
        
        if 0 <= t1 < 1:
            a1 = (1-t1)*(1-t1)*p0s + 2*(1-t1)*t1*p1s + t1*t1*p2s - (r.x*s.x + r.y*s.y)
            a1 /= (s.x*s.x + s.y*s.y)
            if not self.tolerance < a1 <= 1:
                a1 = INFINITY
            if not self.z_height_1 <= a1*dZ <= self.z_height_2:
                a1 = INFINITY
        else:
            a1 = INFINITY
            
        if 0 <= t2 < 1:
            a2 = (1-t2)*(1-t2)*p0s + 2*(1-t2)*t2*p1s + t2*t2*p2s - (r.x*s.x + r.y*s.y)
            a2 /= (s.x*s.x + s.y*s.y)
            if not self.tolerance < a2 <= 1:
                a2 = INFINITY
            if not self.z_height_1 <= a2*dZ <= self.z_height_2:
                a2 = INFINITY
        else:
            a2 = INFINITY
        
        if a2 < a1:
            return a2 * mag_(s)
        else:
            return a1 * mag_(s)
    
    cdef vector_t compute_normal_c(self, vector_t p):
        """use vector calc to get polynomial coefficents then use derivative to find 2d normal
        """
        cdef vector_t normal
        cdef double x1,x2,x3,y1,y2,y3,A,B
        cdef double det,a11,a12,a13,a21,a22,a23
        
        normal.x = 0
        normal.y = 0
        normal.z = 0
        
        x1 = self.X0
        x2 = self.X1
        x3 = self.X2
        y1 = self.Y0
        y2 = self.Y1
        y3 = self.Y2
        
        #solve for A, B and C using these three points,
        #first, find determinant of the three equations (y = Ax^2+Bx+C)
        det = x1*x1*x2 + x1*x3*x3 + x2*x2*x3 - x1*x1*x3 -x1*x2*x2 - x2*x3*x3
        if det == 0:
            print "what the hey? compute normal found det == 0"
            return normal #still just zero's
            
        #then find 6 of the 9 elements of the adjoint matrix (don't have to solve 
        # for C because it is not in derivative
        
        a11 = self.det2x2_(x2,1,x3,1)
        a12 = -self.det2x2_(x1,1,x3,1)
        a13 = self.det2x2_(x1,1,x2,1)
        a21 = -self.det2x2_(x2**2,1,x3**2,1)
        a22 = self.det2x2_(x1**2,1,x3**2,1)
        a23 = -self.det2x2_(x1**2,1,x2**2,1)
        
        #then solve coefficents by multiplying Y by inverse matrix.
        A = y1*a11/det + y2*a12/det +y3*a13/det
        B = y1*a21/det + y2*a22/det +y3*a23/det
        
        normal.y=0      #its a trough shape!
        #reuse variables
        a11 = 2*A*p.x + B
        if a11 == 0:
            normal.z = 1
            normal.x = 0
        else:
            normal.y = -1
            normal.x = a11
        
        #print "debug: calculated bezier normal"
        return norm_(normal)
    
    
cdef int point_in_polygon_c(double X, double Y, object obj):
    cdef int i, size, ct=0
    cdef double y1, y2, h, x, x1, x2
    cdef np.ndarray[np.float64_t, ndim=2] pts=obj
    
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
        if point_in_polygon_c(X,Y, self._xy_points)==1:
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
            
    cdef double intersect_c(self, vector_t p1, vector_t p2):
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
    

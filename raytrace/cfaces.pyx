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
    double INFINITY
    double M_PI
    double sqrt(double)
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)

from ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c

import numpy as np
cimport numpy as np_

cdef struct flatvector_t:
    double x,y


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
    
    
cdef class ElipticalPlaneFace(Face):
    cdef public double g_x, g_y, diameter
    
    params = ['diameter']
    
    def __cinit__(self, **kwds):
        self.g_x = kwds.get('g_x', 0.0)
        self.g_y = kwds.get('g_y', 0.0)
    
    cdef double intersect_c(self, vector_t p1, vector_t p2):
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
        if (X*X + Y*Y) > (d*d/4):
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

#
# Functions used for Bezier math.  should go in utils.pyx when there is one
#


    
cdef double eval_bezier(double t, double cp0, double cp1, double cp2, double cp3):
    #just evaluate a cubic bezier spline
    return cp0*((1-t)**3) + 3*cp1*t*((1-t)**2) + 3*cp2*(1-t)*(t**2) + cp3*(t**3)

cdef double dif_bezier(double t, double cp0, double cp1, double cp2, double cp3):
    #calc the derivative of a cubic bezier when parameter = t
    cdef double A, B, C     #just doin this old school polynomial style
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
     cdef:
      long double a1 = b/a, a2 = c/a, a3 = d/a
      long double Q = (a1*a1 - 3.0*a2)/9.0
      long double R = (2.0*a1*a1*a1 - 9.0*a1*a2 + 27.0*a3)/54.0
      double R2_Q3 = R*R - Q*Q*Q
      double theta
      poly_roots x
     #print "called looking for roots"
     if (R2_Q3 <= 0):
      x.n = 3
      theta = acos(R/sqrt(Q*Q*Q))
      x.roots[0] = -2.0*sqrt(Q)*cos(theta/3.0) - a1/3.0;
      x.roots[1] = -2.0*sqrt(Q)*cos((theta+2.0*M_PI)/3.0) - a1/3.0
      x.roots[2] = -2.0*sqrt(Q)*cos((theta+4.0*M_PI)/3.0) - a1/3.0
     else:
      x.n = 1
      x.roots[0] = pow(sqrt(R2_Q3)+fabs(R), 1/3.0)
      x.roots[0] += Q/x.roots[0]
      x.roots[0] *= (1 if (R < 0.0) else -1)
      x.roots[0] -= a1/3.0
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

        #if p is bigger than atleast one but not all, the it is in the box
        i = p.x > A.x or p.x>B.x or p.x>C.x or p.x>D.x
        j = p.x > A.x and p.x>B.x and p.x>C.x and p.x>D.x
        k = i and not j
        i = p.y > A.y or p.y>B.y or p.y>C.y or p.y>D.y
        j = p.y > A.y and p.y>B.y and p.y>C.y and p.y>D.y
        i = i and not j
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

    cdef double intersect_c(self, vector_t ar, vector_t pee2):

        cdef: 
            flatvector_t tempvector
            flatvector_t r, p2, s, origin
            flatvector_t cp0,cp1,cp2,cp3  #holds control points for spline segment under scrutiny
            double result = INFINITY            #length of ray before it intersects surface. 0 if no valid intersection
            double dZ                       #rate of change of z. dZ*result+Z0 gives Z coordinate
            double A,B,C,D,t,a,b,c,d
            poly_roots ts
            vector_t    tempv
        #print "called intersection"
        origin.x = 0 
        origin.y = 0
        #first off, does ray even enter the depth of the extrusion?
        if (ar.z < self.z_height_1 and pee2.z <self.z_height_1) or (ar.z > self.z_height_2 and pee2.z > self.z_height_2):
            return 0    #does not

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
             ts = roots_of_cubic(A,B,C,D)
             while ts.n > 0:
                 ts.n-=1
                 t = ts.roots[ts.n]
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
                     if b < result and b > self.tolerance:
                      result = b
                
            
        if result == INFINITY: result = 0
        #print "final answer:",result
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
        for curve in self.curves_array.copy():
            cp0.x,cp0.y = curve[0]
            cp1.x,cp1.y = curve[1]
            cp2.x,cp2.y = curve[2]
            cp3.x,cp3.y = curve[3]
            #is point even in this hull?
            if self.pnt_in_hull(ray,cp0,cp1,cp2,cp3):
              #then, solve for t
              #print "in hull"
              cp0 = rotate2D(-theta,cp0)
              cp1 = rotate2D(-theta,cp1)
              cp2 = rotate2D(-theta,cp2)
              cp3 = rotate2D(-theta,cp3)
              #print "cp0: ",cp0.x,cp0.y
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
                if 0<=t<=1:
                 #ok, then is this the t to the same point p? 
                 tmp = eval_bezier(t,cp0.x,cp1.x,cp2.x,cp3.x)
                 #I will generously allow for rounding error 
                 #print "got here with point: ",tmp,
                 if tmp**2 - (ray.x**2+ray.y**2) < .005:
                    #this is the single solution. return the derivative dy/dx = dy/dt / dx/dt
                    
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
        print "error: Bezier normal not found, point not actually on curve!"
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
    

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

"""
Cython module for Face definitions
"""
cdef extern from "math.h":
    double INFINITY

from ctracer cimport Face, sep_, intersection_t,\
        vector_t, ray_t, FaceList


cdef class CircularFace(Face):
    cdef public double diameter, offset, z_plane
    
    params = ['diameter', 'offset']
    
    def __cinit__(self, **kwds):
        self.z_plane = kwds.get('z_plane', 0.0)
    
    cdef intersection_t intersect_c(self, vector_t p1, vector_t p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        cdef:
            double max_length, h, X, Y
            intersection_t inter
            vector_t point
        #print "CFACE", p1, p2
        max_length = sep_(p1, p2)
        h = (self.z_plane-p1.z)/(p2.z-p1.z)
        if (h<self.tolerance) or (h>1.0):
            #print "H", h
            inter.dist = INFINITY
            return inter
        X = p1.x + h*(p2.x-p1.x) - self.offset
        Y = p1.y + h*(p2.y-p1.y)
        if (X**2 + Y**2) > (self.diameter/2.)**2:
            #print "X", X, "Y", Y
            inter.dist = INFINITY
            return inter
        inter.dist = max_length * h
        point.x = X + self.offset
        point.y = Y
        point.z = self.z_plane
        inter.point = point
        inter.face_idx = self.idx
        #print "CF_INTER", inter.point, inter.dist, inter.face_idx
        return inter

    cdef vector_t compute_normal_c(self, vector_t p):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef vector_t normal
        normal.x=0
        normal.y=0
        normal.z=-1
        return normal
    
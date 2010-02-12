"""
Cython module for Face definitions
"""
cdef extern from "math.h":
    double INFINITY

from ctracer cimport Face, sep_, eval_PEC_children, intersection_t,\
        vector_t, ray_t, FaceList


cdef class CircularFace(Face):
    cdef public double diameter, offset
    
    params = ['diameter', 'offset']
    
    cdef intersection_t intersect_c(self, vector_t p1, vector_t p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        cdef:
            double max_length, h, X, Y
            intersection_t inter
            vector_t point
            
        max_length = sep_(p1, p2)
        h = -p1.z/(p2.z-p1.z)
        if (h<self.tolerance) or (h>1.0):
            inter.dist = INFINITY
            return inter
        X = p1.x + h*(p2.x-p1.x) - self.offset
        Y = p1.y + h*(p2.y-p1.y)
        if (X**2 + Y**2) > (self.diameter/2.)**2:
            inter.dist = INFINITY
            return inter
        inter.dist = max_length * h
        point.x = X + self.offset
        point.y = Y
        point.z = 0.0
        inter.point = point
        inter.face_idx = self.idx
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
    
    cdef ray_t eval_child_ray_c(self, ray_t old_ray, int ray_idx, 
                                vector_t p, FaceList face_set):
        return eval_PEC_children(<Face>self, old_ray, ray_idx,
                                p, face_set)
"""
Cython module for Face definitions
"""

from ctracer cimport Face, sep_

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
    
#    def intersect(self, P1, P2, max_lenth):
#        """
#        Calculate intersection point for a ray segment between two points
#        P1 and P2, in the local coordinate system.
#        """
#        max_length = numpy.sqrt(((P2 - P1)**2).sum(axis=1))
#        r = self.diameter/2
#        offset = self.offset
        
#        x1,y1,z1 = P1.T
#        x2,y2,z2 = P2.T
        
#        h = -z1/(z2-z1)
#        X = x1 + h*(x2-x1) - offset
#        Y = y1 + h*(y2-y1)
        
#        length = max_length*h
        
#        mask = (X**2 + Y**2) > r**2
#        mask = numpy.logical_or(mask, h < self.tolerance)
#        mask = numpy.logical_or(mask, h > 1.0)
        
#        length[mask] = numpy.Infinity
        
#        t_points = numpy.column_stack((X+offset, Y, numpy.zeros_like(X)))
    
#        dtype=([('length','f8'),('face', 'O'),('point','f8',3)])
#        result = numpy.empty(P1.shape[0], dtype=dtype)
#        result['length'] = length
#        result['face'] = self
#        result['point'] = t_points
#        return result

    cdef vector_t compute_normal_c(self, vector_t p):
        cdef vector_t normal
        return rotate_c(self.trans
    
#    def compute_normal(self, points):
#        """computes the surface normal in the global coordinate system"""
#        t = self.transform
#        n = numpy.asarray([t.transform_vector(0,0,-1),])
#        return numpy.ones(points.shape) * n
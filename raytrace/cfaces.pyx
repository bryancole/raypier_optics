"""
Cython module for Face definitions
"""
cdef extern from "listobject.h":
    ctypedef class __builtin__.list [object PyListObject]:
        pass


#cdef class Matrix(object):

cdef struct t_vector:
    float x
    float y
    float z
    
cdef struct Complex:
    float real
    float imag
    
    
cdef class Ray:
    cdef t_vector origin, direction, E_vector
    cdef public float length, offset_length, E1_amp, E2_amp
    cdef Complex refractive_index
    
    def __cinit__(self, origin=[0.0,0.0,0.0], 
                    direction=[0.0,0.0,0.0], 
                    e_vec=[0.0,0.0,0.0]):
        self.origin.x = origin[0]
        self.origin.y = origin[1]
        self.origin.z = origin[2]
        self.direction.x = direction[0]
        self.direction.y = direction[1]
        self.direction.z = direction[2]
        self.e_vec.x = e_vec[0]
        self.e_vec.y = e_vec[1]
        self.e_vec.z = e_vec[2]
    
    
cdef class Vector(object):
    cdef public float x, y, z
    
    def __cinit__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z
    
    def __repr__(self):
        return "%s(%g, %g, %g)"%(self.__class__.__name__,
                                self.x, self.y, self.z)
    
    
cdef class Point(Vector):
    pass
    

cdef class Face(object):
    cdef public object owner
    cdef public char *name
    cdef public float tolerance
    
    def __cinit__(self):
        self.name = "base Face class"
        self.tolerance = 0.0001
    
    cpdef float intersect(self, Point p1, Point p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        return 0.5

    cpdef Vector compute_normal(self, Point p):
        return Vector(p.z,p.y,p.x)
    
    cpdef Ray eval_child_ray(self, Ray ray, Point p):
        return ray
    
    
cdef class Intersection:
    cdef Face face
    cdef Point point
    
    def __cinit__(self, face, point):
        self.face = face
        self.point = point
    
    
cdef class FaceList(list):
    """A group of faces which share a transform"""
    cdef public object transform
    
    cpdef Vector transform(self, Vector v):
        return v
    
#    cpdef Vector inv_transform(self, Vector v):
#        return v
#    
#    cpdef Intersection intersect(self, Point P1, Point P2):
#        """Finds the face with the nearest intersection
#        point, for the ray defined by the two input points,
#        P1 and P2 (in global coords).
#        """
#        return Intersection(self[0], P1)

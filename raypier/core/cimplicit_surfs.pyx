from scipy.ndimage._ni_docstrings import _origin_doc

cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)
    
from .ctracer cimport ImplicitSurface, sep_, \
        vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, cross_
            
import numpy as np
cimport numpy as np_

from .ctracer import ImplicitSurface


cdef class NullSurface(ImplicitSurface):
    cdef double evaluate_c(self, vector_t p) nogil:
        return -1.0


cdef class Plane(ImplicitSurface):
    cdef:
        vector_t _origin
        vector_t _normal
        
    def __cinit__(self, **kwds):
        self.origin = kwds.get('origin', (0.0,0.0,0.0))
        self.normal = kwds.get('normal', (1.0,0.0,0.0))
        
    property origin:
        def __set__(self, o):
            self._origin.x = o[0]
            self._origin.y = o[1]
            self._origin.z = o[2]
            
        def __get__(self):
            cdef vector_t o=self._origin
            return (o.x, o.y, o.z)
        
    property normal:
        def __set__(self, n):
            cdef:
                vector_t _n
            _n.x = n[0]
            _n.y = n[1]
            _n.z = n[2]
            self._normal = norm_(_n)
            
        def __get__(self):
            cdef vector_t n=self._normal
            return (n.x, n.y, n.z)
        
    cdef double evaluate_c(self, vector_t p) nogil:
        return dotprod_(self._normal, subvv_(p, self._origin))
    

cdef class Sphere(ImplicitSurface):
    cdef:
        vector_t _centre
        public double radius
        
    def __cinit__(self, **kwds):
        self.centre = kwds.get('centre', (0.0,0.0,0.0))
        self.radius = kwds.get('radius', 1.0)
        
    property centre:
        def __set__(self, o):
            self._centre.x = o[0]
            self._centre.y = o[1]
            self._centre.z = o[2]
            
        def __get__(self):
            cdef vector_t o=self._centre
            return (o.x, o.y, o.z)
        
    cdef double evaluate_c(self, vector_t p) nogil:
        return sep_(p, self._centre) - self.radius
    
    
cdef class Cylinder(ImplicitSurface):
    cdef:
        vector_t _origin
        vector_t _axis
        public double radius
        
    def __cinit__(self, **kwds):
        self.origin = kwds.get('origin', (0.0,0.0,0.0))
        self.axis = kwds.get('axis', (0.,0.,1.))
        self.radius = kwds.get('radius', 10.0)
        
    property origin:
        def __get__(self):
            cdef vector_t o=self._origin
            return (o.x, o.y, o.z)
        
        def __set__(self, v):
            cdef vector_t o
            o.x = v[0]
            o.y = v[1]
            o.z = v[2]
            self._origin = o
            
    property axis:
        def __get__(self):
            cdef vector_t o=self._axis
            return (o.x, o.y, o.z)
        
        def __set__(self, v):
            cdef vector_t o
            o.x = v[0]
            o.y = v[1]
            o.z = v[2]
            self._axis = norm_(o)
            
    cdef double evaluate_c(self, vector_t p) nogil:
        return mag_(cross_(subvv_(p, self._origin), self._axis)) - self.radius
    

cdef class Union(ImplicitSurface):
    cdef:
        public np_.ndarray _surfaces
        void** _surfaces_ptr
        int size
        
    def __cinit__(self, *args):
        self.surfaces = args
        
    property surfaces:
        def __get__(self):
            return self._surfaces
        
        def __set__(self, oa):
            self._surfaces = np.asarray(oa, dtype='O')
            self.size = self._surfaces.shape[0]
            self._surfaces_ptr = <void**>(self._surfaces.data)
            
        
#     property surfaces:
#         def __set__(self, surf_list):
#             if not all(isinstance(s,ImplicitSurface) for s in surf_list):
#                 raise ValueError("All surfaces must be instances of ImplicitSurface")
#             self._surfaces = list(surf_list)
#             
#         def __get__(self):
#             return self._surfaces
        
    cdef double evaluate_c(self, vector_t p) nogil:
        cdef:
            void* surf
            int i
            double out, v
            
        surf = self._surfaces_ptr[0]
        out = (<ImplicitSurface>surf).evaluate_c(p)
        for i in range(1, self.size):
            surf = self._surfaces_ptr[i]
            v = (<ImplicitSurface>surf).evaluate_c(p)
            out = self.apply_op(out, v)
        return out
    
    cdef double apply_op(self, double out, double v) nogil:
        if v < out:
            out = v
        return out
    
    
cdef class Intersection(Union):
    cdef double apply_op(self, double out, double v) nogil:
        if v > out:
            out = v
        return out
    

cdef class Difference(Union):
    cdef double apply_op(self, double out, double v) nogil:
        out -= v
        return out

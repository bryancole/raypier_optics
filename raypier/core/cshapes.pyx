
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)


from .ctracer cimport Shape, sep_, \
        vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, cross_
            
import numpy as np
            

cdef class LogicalOpShape(Shape):
    def __and__(self, Shape other):
        return BooleanAND(self, other)
    
    def __or__(self, Shape other):
        return BooleanOR(self, other)
    
    def __xor__(self, Shape other):
        return BooleanXOR(self, other)
    
    def __invert__(self):
        return InvertShape(self)
            
            
cdef class InvertShape(LogicalOpShape):
    cdef:
        public Shape shape
        
    def __cinit__(self, Shape shape):
        self.shape = shape
        
    cdef bint point_inside_c(self, double x, double y):
        return 1 & (~self.shape.point_inside_c(x,y))
            
            
cdef class BooleanShape(LogicalOpShape):
    cdef:
        public Shape shape1
        public Shape shape2
        
    def __cinit__(self, Shape shape1, Shape shape2):
        self.shape1 = shape1
        self.shape2 = shape2
        
        
cdef class BooleanAND(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return (<Shape>self.shape1).point_inside_c(x,y) & (<Shape>self.shape2).point_inside_c(x,y)
    
    
cdef class BooleanOR(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return (<Shape>self.shape1).point_inside_c(x,y) | (<Shape>self.shape2).point_inside_c(x,y)
    
    
cdef class BooleanXOR(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return (<Shape>self.shape1).point_inside_c(x,y) ^ (<Shape>self.shape2).point_inside_c(x,y)
    
    
cdef class BasicShape(LogicalOpShape):
    cdef:
        double centre_x, centre_y
        
    def __cinit__(self, **kwds):
        if "centre" in kwds:
            self.centre = kwds["centre"]
        else:
            self.centre_x = kwds.get("centre_x", 0.0)
            self.centre_y = kwds.get("centre_y", 0.0)
        
    property centre:
        def __get__(self):
            return (self.centre_x, self.centre_y)
        
        def __set__(self, v):
            self.centre_x = v[0]
            self.centre_y = v[1]
            
    def __and__(self, Shape other):
        return BooleanAND(self, other)
    
    def __or__(self, Shape other):
        return BooleanOR(self, other)
    
    def __xor__(self, Shape other):
        return BooleanXOR(self, other)
    
    def __invert__(self):
        return InvertShape(self)
                
            
cdef class CircleShape(BasicShape):
    cdef:
        public double radius
        
    def __cinit__(self, **kwds):
        self.radius = kwds.get("radius", 1.0)
        
    cdef bint point_inside_c(self, double x, double y):
        cdef:
            double dx = x-self.centre_x
            double dy = y-self.centre_y
        if (dx*dx) + (dy*dy) < (self.radius*self.radius):
            return 1
        else:
            return 0
        
        
cdef class RectangleShape(BasicShape):
    cdef:
        public double width
        public double height
        
    def __cinit__(self, **kwds):
        self.width = kwds.get("width", 5.0)
        self.height = kwds.get("height", 7.0)
        
    cdef bint point_inside_c(self, double x, double y):
        cdef:
            double dx = x-self.centre_x
            double dy = y-self.centre_y
            
        if (2*fabs(dx) < self.width) and (2*fabs(dy)) < self.height:
            return 1
        else:
            return 0


cdef class PolygonShape(BasicShape):
    cdef:
        double[:,:] _coordinates
        
    property coordinates:
        def __get__(self):
            return np.asarray(self._coordinates)
        
        def __set__(self, val):
            self._coordinates = val
            
    cdef bint point_inside_c(self, double X, double Y):
        cdef:
            int i, size, ct=0
            double y1, y2, x1, x2
            double[:,:] pts = self._coordinates
        
        size = pts.shape[0]
        
        y1 = pts[size-1,1]
        x1 = pts[size-1,0]
        for i in range(size):
            y2 = pts[i,1]
            x2 = pts[i,0]
            if (y1 <= Y < y2) or (y2 <= Y < y1):
                if (x1 + (Y-y1)*(x2 - x1)/(y2 - y1)) > X:
                    ct = not ct
            y1 = y2
            x1 = x2
        return ct
    
        
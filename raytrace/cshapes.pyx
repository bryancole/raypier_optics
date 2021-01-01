
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)


from ctracer cimport Shape, sep_, \
        vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, cross_
            
            
cdef class InvertShape(Shape):
    cdef:
        public Shape shape
        
    def __cinit__(self, Shape shape):
        self.shape = shape
        
    cdef bint point_inside_c(self, double x, double y):
        return 1 & (~self.shape.point_inside_c(x,y))
            
            
cdef class BooleanShape(Shape):
    cdef:
        public Shape shape1
        public Shape shape2
        
    def __cinit__(self, Shape shape1, Shape shape2):
        self.shape1 = shape1
        self.shape2 = shape2
        
        
cdef class BooleanAND(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return self.shape1.point_inside_c(x,y) & self.shape2.point_inside_c(x,y)
    
    
cdef class BooleanOR(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return self.shape1.point_inside_c(x,y) | self.shape2.point_inside_c(x,y)
    
    
cdef class BooleanXOR(BooleanShape):
    cdef bint point_inside_c(self, double x, double y):
        return self.shape1.point_inside_c(x,y) ^ self.shape2.point_inside_c(x,y)
    
    
cdef class BasicShape(Shape):
    cdef:
        double centre_x, centre_y
        
    def __cinit__(self, **kwds):
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
        self.centre_x = kwds.get("radius", 1.0)
        
    cdef bint point_inside_c(self, double x, double y):
        cdef:
            double dx = x-self.centre_x
            double dy = y-self.centre_y
        if (dx*dx) + (dy*dy) < (self.radius*self.radius):
            return 1
        else:
            return 0
        
        

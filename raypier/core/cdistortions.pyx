
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x ) nogil
    double pow(double x, double y) nogil
    double fabs(double) nogil
    double cos(double) nogil
    double sin(double) nogil
    double acos(double) nogil

cdef extern from "float.h":
    double DBL_MAX

cdef double INF=(DBL_MAX+DBL_MAX)

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, Shape, Distortion
                
cimport numpy as np_
import numpy as np
                
                
cdef class SimpleTestZernikeJ7(Distortion):
    """Implement one low-order Zernike poly for testing purposes.
    """
    cdef:
        public double unit_radius
        public double amplitude
        
    def __cinit__(self, **kwds):
        self.unit_radius = kwds.get("unit_radius", 1.0)
        self.amplitude = kwds.get("amplitude", 1.0)
        
    cdef double z_offset_c(self, double x, double y) nogil:
#         cdef:
#             double rho, Z
            
        x /= self.unit_radius
        y /= self.unit_radius
        #rho = sqrt(x*x + y*y)
        #theta = atan2(y, x)
        Z = sqrt(8.0)*(3*(x*x + y*y) - 2)*y
        return Z * self.amplitude
    
    cdef vector_t z_offset_and_gradient_c(self, double x, double y) nogil:
        """The z-axis surface sag is returned as the z-component 
        of the output vector. The x- and y-components of the surface
        gradient are placed in the x- and y- components of vector.
        I.e. the returned vector is
            (dz(x,y)/dx, dz(x,y)/dy, z(x,y))
        """
        cdef:
            double rho, theta, Z, root8=sqrt(8)*self.amplitude, R=self.unit_radius
            vector_t p
            
        x /= R
        y /= R
            
        #rho = sqrt(x*x + y*y)
        #theta = atan2(y, x)
            
        p.x = root8 * 6*x*y /R #root8*(6*y + (3*rho*rho-2)/y)*x
        p.y = root8 * (3*x*x + 9*y*y -2)/R#root8*(6*y*y + 3*rho*rho - 2)
        p.z = root8*(3*(x*x + y*y) - 2)*y 
        return p
        
                
cdef class ZernikeDistortion(Distortion):
    cdef:
        double[:] _coefs
         
    def __cinit__(self, coefs):
        self._coefs = np.ascontiguousarray(coefs, 'd')
         
        nk=len(coefs)
        n = int(np.ceil((np.sqrt(8*nk+1)-3)/2))
        nk = (n+1)*(n+2)//2
        
        
        
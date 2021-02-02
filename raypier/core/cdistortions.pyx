
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)

cdef extern from "float.h":
    double DBL_MAX

cdef double INF=(DBL_MAX+DBL_MAX)

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, Shape, Distortion
                
cimport numpy as np_
import numpy as np
                
                
                
cdef class ZernikeDistortion(Distortion):
    cdef:
        double[:] _coefs
        
    def __cinit__(self, coefs):
        self._coefs = np.ascontiguousarray(coefs, 'd')
        
        nk=len(coefs)
        n = int(np.ceil((np.sqrt(8*nk+1)-3)/2))
        nk = (n+1)*(n+2)//2
        
        
        
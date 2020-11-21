from IPython.testing import iptest

cdef extern from "math.h":
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

from ctracer import ray_dtype
import numpy as np
cimport numpy as np_


cpdef  np_.npy_complex128[:,:] sum_gaussian_modes(ray_t[:] rays, np_.npy_float64[:] phases, 
                       np_.npy_float64[:,:] modes, 
                       np_.npy_float64[:,:] points):
    """
    Compute the E-field at the given points by propagation of the
    given rays, with phase and mode-coefficients.
    """
    cdef:
        size_t iray, ipt 
        size_t Nray=rays.shape[0]
        size_t Npt=points.shape[0]
        
        ray_t ray
        np_.npy_complex128[:,:] out = np.empty((Npt,3), dtype=np.complex128)
        
    for iray in range(Nray):
        ray = rays[iray]
        for ipt in range(Npt):
            out[ipt,0].real += ray.origin.x
            out[ipt,0].imag += ray.direction.x
            out[ipt,1].real += ray.origin.y
            out[ipt,1].imag += ray.direction.y
            out[ipt,2].real += ray.origin.z
            out[ipt,2].imag += ray.direction.z
            
    return out

        
    
    
    
    
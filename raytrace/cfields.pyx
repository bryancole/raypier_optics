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
    

IF UNAME_SYSNAME == "Windows":
    cdef extern from "complex.h":
        double complex csqrt "sqrt" (double complex)
        double cabs "abs" (double complex)
        double complex cexp "exp" (double complex)
    
    cdef double complex I = 1j
ELSE:
    cdef extern from "complex.h":
        double complex csqrt (double complex)
        double cabs (double complex)
        double complex cexp (double complex)
        double complex I    


from ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, cross_

from ctracer import ray_dtype
import numpy as np
cimport numpy as np_


cpdef  np_.npy_complex128[:,:] sum_gaussian_modes(ray_t[:] rays, np_.npy_float64[:] phases, 
                       double complex[:,:] modes, np_.npy_float64[:] wavelengths,
                       vector_t[:] points):
    """
    Compute the E-field at the given points by propagation of the
    given rays, with phase and mode-coefficients.
    """
    cdef:
        size_t iray, ipt 
        size_t Nray=rays.shape[0]
        size_t Npt=points.shape[0]
        
        ray_t ray
        np_.npy_complex128[:,:] out = np.zeros((Npt,3), dtype=np.complex128)
        vector_t pt, H, E
        double complex U, A, B, C, detG0, denom, AA, CC
        double x,y,z, phase, k
        
    for iray in range(Nray):
        
        ray = rays[iray]
        E = ray.E_vector
        H = cross_(E, ray.direction)
        A = modes[iray, 0]
        B = modes[iray, 1]
        C = modes[iray, 2]
        detG0 = (A*C) - (B*B)
        phase = phases[iray]
        k = 2*M_PI/wavelengths[ray.wavelength_idx]
        
        for ipt in range(Npt):
            pt = subvv_(points[ipt], ray.origin)
            x = dotprod_(pt, E)
            y = dotprod_(pt, H)
            z = dotprod_(pt, ray.direction)
            
            denom = 1 + (z*(A+C)) + (z**2)*detG0
            AA = (A + z*detG0)/denom
            CC = (C + z*detG0)/denom
            U = cexp((I*phase) + I*k*(z +  AA*(x**2) + (2*B*x*y) + CC*(y**2) ))
            ###Normalisation factor
            U /= csqrt((1 + z*A)*(1 + z*C) - (z*B)*(z*B))
            
            ###reuse AA and CC
            AA.real = ray.E1_amp.real
            AA.imag = ray.E1_amp.imag
            AA *= U
            CC.real = ray.E2_amp.real
            CC.imag = ray.E2_amp.imag
            CC *= U
            
            out[ipt,0].real += (AA.real * E.x) + (CC.real * H.x) 
            out[ipt,0].imag += (AA.imag * E.x) + (CC.imag * H.x)
            out[ipt,1].real += (AA.real * E.y) + (CC.real * H.y)
            out[ipt,1].imag += (AA.imag * E.y) + (CC.imag * H.y)
            out[ipt,2].real += (AA.real * E.z) + (CC.real * H.z)
            out[ipt,2].imag += (AA.imag * E.z) + (CC.imag * H.z)
            
    return out

        
    
    
    
    
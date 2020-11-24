
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    double cos(double)
    double sin(double)
    double acos(double)
    

IF UNAME_SYSNAME == "Windows":
    cdef extern from "complex.h":
        double complex csqrt "sqrt" (double complex) nogil
        double cabs "abs" (double complex) nogil
        double complex cexp "exp" (double complex) nogil
    
    cdef double complex I = 1j
ELSE:
    cdef extern from "complex.h":
        double complex csqrt (double complex) nogil
        double cabs (double complex) nogil
        double complex cexp (double complex) nogil
        double complex I    
        
from cython.parallel import prange


from ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, cross_, RayCollection

from ctracer import ray_dtype
cimport cython
import numpy as np
cimport numpy as np_


cpdef check_this(ray_t[:] rays):
    return rays.shape


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.cdivision(True)
cpdef  sum_gaussian_modes(RayCollection rays,
                          double complex[:,:] modes, 
                          np_.npy_float64[:] wavelengths,
                          np_.npy_float64[:,:] points):
    """
    Compute the E-field at the given points by propagation of the
    given rays, with phase and mode-coefficients.
    """
    cdef:
        size_t iray, ipt 
        size_t Nray=rays.n_rays
        size_t Npt=points.shape[0]
        
        ray_t ray
        np_.npy_complex128[:,:] out = np.zeros((Npt,3), dtype=np.complex128)
        vector_t pt, H, E
        double complex U, A, B, C, E1, E2, detG0
        double x,y,z, phase, k, inv_root_area
    
    with nogil:
        for iray in range(Nray):
            
            ray = rays.rays[iray]
            E = ray.E_vector
            H = cross_(E, ray.direction)
            k = 2000.0*M_PI/wavelengths[ray.wavelength_idx]
            A = modes[iray, 0]
            B = modes[iray, 1]
            C = modes[iray, 2]
            A.imag /= k
            B.imag /= k
            C.imag /= k
            detG0 = (A*C) - (B*B)
            phase = ray.phase
            
            ###normalisation factor 1/root(area)
            inv_root_area = sqrt(sqrt(A.imag*C.imag -(B.imag*B.imag)))/M_PI
            
            for ipt in prange(Npt):
                pt.x = points[ipt,0]
                pt.y = points[ipt,1]
                pt.z = points[ipt,2]
                pt = subvv_(pt, ray.origin)
                
                U = calc_mode_U(A,B,C, detG0, pt, E, H, ray.direction, k, phase, inv_root_area)
                
                E1 = (ray.E1_amp.real + I*ray.E1_amp.imag) * U
                E2 = (ray.E2_amp.real + I*ray.E2_amp.imag) * U
                 
                out[ipt,0].real += (E1.real * E.x) + (E2.real * H.x) 
                out[ipt,0].imag += (E1.imag * E.x) + (E2.imag * H.x)
                out[ipt,1].real += (E1.real * E.y) + (E2.real * H.y)
                out[ipt,1].imag += (E1.imag * E.y) + (E2.imag * H.y)
                out[ipt,2].real += (E1.real * E.z) + (E2.real * H.z)
                out[ipt,2].imag += (E1.imag * E.z) + (E2.imag * H.z)
            
    return np.asarray(out)

        
cdef double complex calc_mode_U(double complex A,
                                double complex B,
                                double complex C,
                                double complex detG0,
                                vector_t pt, 
                                vector_t E, 
                                vector_t H,
                                vector_t direction,
                                double k,
                                double phase, 
                                double inv_root_area) nogil:
    cdef:
        double complex denom, AA, CC, U
        double x,y,z
        
    x = dotprod_(pt, E)
    y = dotprod_(pt, H)
    z = dotprod_(pt, direction)
    denom = 1 + (z*(A+C)) + (z**2)*detG0
    AA = (A + z*detG0)/denom
    CC = (C + z*detG0)/denom
    U = cexp( (I*phase) + I*k*(z + AA*(x**2) + (2*B*x*y)/denom + CC*(y**2) ) )
    ###Normalisation factor
    U /= csqrt((1 + z*A)*(1 + z*C) - (z*B)*(z*B))
    
    ###normalise by ray initial area
    U *= inv_root_area
    
    return U
    
    
    
    
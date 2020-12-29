
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

cdef:
    double complex rootI=csqrt(I)


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
            E = norm_(ray.E_vector)
            H = norm_(cross_(ray.direction, E))
            k = 2000.0*M_PI/wavelengths[ray.wavelength_idx]
            A = modes[iray, 0]
            B = modes[iray, 1]
            C = modes[iray, 2]
            
            ###normalisation factor 1/root(area)
            inv_root_area = sqrt(sqrt(A.imag*C.imag -(B.imag*B.imag))/M_PI)
            
            A.imag /= k
            B.imag /= k
            C.imag /= k
            detG0 = (A*C) - (B*B)
            phase = ray.phase + (ray.accumulated_path*k)
            
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
    denom = 1 + (z*(A+C)) + (z*z)*detG0
    AA = (A + z*detG0)/denom
    CC = (C + z*detG0)/denom
    U = cexp( (I*phase) + I*k*(z + AA*(x*x) + B*(2*x*y)/denom + CC*(y*y) ) )
    ###Normalisation factor
    ### I've added an extra 1I factor here in order to rotate the phase away from the 0/2pi wrapping point.
    ### Seems like the csqrt (or cpow) functions struggle with accuracy near phase zero where we're summing 
    ### many values straddling the wrapping point
    U /= csqrt(((1 + z*A)*(1 + z*C) - (z*B)*(z*B))*I)
    
    ###normalise by ray initial area
    ### multiply by sqrt(I) to reverse the effect of the I factor applied above
    U *= inv_root_area*rootI
    
    return U


cdef void evaluate_one_mode(np_.npy_complex128[:] out, double[:] x, double[:] y, double[:] dx, double[:] dy, double blending, int size):
    ### Brute force evaluation of the linear least squares fit to the x,y,dx,dy coefs
    ### We solve A.T*A x = A.T b which (by magical means) gives the least squares fit to A x = b
    ### A.T*A and A.T b have 3 rows i.e. we just need to solve 3 equations for 3 unknowns.
    ### Even better, the matrices are symmetric. Solving in sympy gives the closed form of the solution
    ### as fairly simple expressions 
    cdef:
        double a00=0, a01=0, a02=0, a11=0, a12=0, a22=0, b0=0, b1=0, b2=0
        int i
        double xi2, yi2
    
    ### Imag parts ###    
    for i in range(size):
        xi2 = x[i]*x[i]
        yi2 = y[i]*y[i]
        a00 += xi2*xi2
        a01 += 2*xi2*x[i]*y[i]
        a02 += xi2*yi2
        a11 += 4*xi2*yi2
        a12 += 2*x[i]*y[i]*yi2
        a22 += yi2*yi2
        b0 += xi2
        b1 += 2*x[i]*y[i]
        b2 += yi2
    
    out[0].imag = blending *(-b0*(a11*a22 - a12**2) + b1*(a01*a22 - a02*a12) - b2*(a01*a12 - a02*a11))/(-a00*a11*a22 + a00*a12**2 + a01**2*a22 - 2*a01*a02*a12 + a02**2*a11)
    out[1].imag = blending *(b0*(a01*a22 - a02*a12) - b1*(a00*a22 - a02**2) + b2*(a00*a12 - a01*a02))/(-a00*a11*a22 + a00*a12**2 + a01**2*a22 - 2*a01*a02*a12 + a02**2*a11)
    out[2].imag = blending *(-b0*(a01*a12 - a02*a11) + b1*(a00*a12 - a01*a02) - b2*(a00*a11 - a01**2))/(-a00*a11*a22 + a00*a12**2 + a01**2*a22 - 2*a01*a02*a12 + a02**2*a11)
    
    ### Real parts ###
    a00=0
    a01=0
    a02=0
    a11=0
    a12=0
    a22=0
    b0=0
    b1=0
    b2=0 
    for i in range(size):
        xi2 = x[i]*x[i]
        yi2 = y[i]*y[i]
        
        a00 += xi2
        a01 += x[i]*y[i]
        a11 += xi2 + yi2
        a22 += yi2
        
        b0 += dx[i]*x[i]
        b1 += dx[i]*y[i] + dy[i]*x[i]
        b2 += dy[i]*y[i]
        
    a12 = a01*a01 #reuse a12 as an optimisation

    out[0].real = (-a12*b2 + a01*a22*b1 + b0*(a12 - a11*a22))/(a00*a12 - a00*a11*a22 + a12*a22)
    out[1].real = (a00*a01*b2 - a00*a22*b1 + a01*a22*b0)/(a00*a12 - a00*a11*a22 + a12*a22)
    out[2].real = (a00*a01*b1 - a12*b0 - b2*(a00*a11 - a12))/(a00*a12 - a00*a11*a22 + a12*a22)
    
    
def evaluate_modes(double[:,:] x, double[:,:] y, double[:,:] dx, double[:,:] dy, double blending=1.0):
    cdef:
        int n_rays=x.shape[0], row_size=x.shape[1]
        int i
        np_.npy_complex128[:,:] out = np.zeros((n_rays,3), dtype=np.complex128)        
    
    for i in range(n_rays):
        evaluate_one_mode(out[i,:], x[i], y[i], dx[i], dy[i], blending, row_size)
        
    return out
    

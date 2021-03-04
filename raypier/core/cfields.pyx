
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


from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, cross_, RayCollection, complex_t, \
                GaussletCollection, gausslet_t, para_t

from .ctracer import ray_dtype
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
        complex_t[:,:] out = np.zeros((Npt,3), dtype=np.complex128)
        vector_t pt, H, E
        double complex U, A, B, C, E1, E2, detG0, kz
        double x,y,z, phase, k, inv_root_area, invk
    
    with nogil:
        for iray in range(Nray):
            ray = rays.rays[iray]
            E = norm_(ray.E_vector)
            H = norm_(cross_(ray.direction, E))
            k = 2000.0*M_PI/wavelengths[ray.wavelength_idx]
            ### The accumulated path length already includes the refractive index of up-stream rays
            ### Hence, need to calculate it before applying the refractive index of the last leg.
            phase = ray.phase + (ray.accumulated_path*k)
            
            #kz = ray.refractive_index.real + I*ray.refractive_index.imag
            kz = ray.refractive_index
            kz *= k
            invk = 2./kz.real
            A = modes[iray, 0]
            B = modes[iray, 1]
            C = modes[iray, 2]
            
            ###normalisation factor 1/root(area). Yes, there really is a double square root.
            inv_root_area = sqrt(sqrt(A.imag*C.imag -(B.imag*B.imag))*(2.0/M_PI))
            
            A.imag *= invk
            B.imag *= invk
            C.imag *= invk
            detG0 = (A*C) - (B*B)
            
            for ipt in prange(Npt):
                pt.x = points[ipt,0]
                pt.y = points[ipt,1]
                pt.z = points[ipt,2]
                pt = subvv_(pt, ray.origin)
                
                U = calc_mode_U(A,B,C, detG0, pt, E, H, ray.direction, kz, phase, inv_root_area)
                
                E1 = ray.E1_amp * U
                E2 = ray.E2_amp * U
                
                out[ipt,0] += (E1*E.x + E2*H.x)
                out[ipt,1] += (E1*E.y + E2*H.y)
                out[ipt,2] += (E1*E.z + E2*H.z)
    return np.asarray(out)

        
cdef double complex calc_mode_U(double complex A,
                                double complex B,
                                double complex C,
                                double complex detG0,
                                vector_t pt, 
                                vector_t E, 
                                vector_t H,
                                vector_t direction,
                                double complex k,
                                double phase, 
                                double inv_root_area) nogil:
    cdef:
        double complex denom, AA, CC, U
        double x,y,z
        
    x = dotprod_(pt, E)
    y = dotprod_(pt, H)
    z = dotprod_(pt, direction)
    ### Where did this factor of 2 come from?
    denom = (1 + (z*(A+C)) + (z*z)*detG0) *2
    AA = (A + z*detG0)/denom
    CC = (C + z*detG0)/denom
    U = cexp( I*(phase + k*(z + AA*(x*x) + B*(2*x*y)/denom + CC*(y*y) ) ) )
    ###Normalisation factor
    ### I've added an extra 1I factor here in order to rotate the phase away from the 0/2pi wrapping point.
    ### Seems like the csqrt (or cpow) functions struggle with accuracy near phase zero where we're summing 
    ### many values straddling the wrapping point
    U /= csqrt(((1 + z*A)*(1 + z*C) - (z*B)*(z*B))*I)
    
    ###normalise by ray initial area
    ### multiply by sqrt(I) to reverse the effect of the I factor applied above
    U *= inv_root_area*rootI
    
    return U


@cython.boundscheck(False)
@cython.wraparound(False)
cdef void evaluate_one_mode(np_.npy_complex128[:] out, double[:] x, double[:] y, 
                            double[:] dx, double[:] dy, double blending, int size) nogil:
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
    out[1].imag = -blending *(b0*(a01*a22 - a02*a12) - b1*(a00*a22 - a02**2) + b2*(a00*a12 - a01*a02))/(-a00*a11*a22 + a00*a12**2 + a01**2*a22 - 2*a01*a02*a12 + a02**2*a11)
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
    out[1].real = -(a00*a01*b2 - a00*a22*b1 + a01*a22*b0)/(a00*a12 - a00*a11*a22 + a12*a22)
    out[2].real = (a00*a01*b1 - a12*b0 - b2*(a00*a11 - a12))/(a00*a12 - a00*a11*a22 + a12*a22)
    
    return
    
    
def evaluate_modes(double[:,:] x, double[:,:] y, double[:,:] dx, double[:,:] dy, double blending=1.0):
    cdef:
        int n_rays=x.shape[0], row_size=x.shape[1]
        int i
        np_.npy_complex128[:,:] out = np.zeros((n_rays,3), dtype=np.complex128)        
    
    with nogil:
        for i in prange(n_rays):
            evaluate_one_mode(out[i,:], x[i], y[i], dx[i], dy[i], blending, row_size)
        
    return np.asarray(out)

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def calc_mode_curvature(double[:] rx, double[:] ry, double[:] dx, double[:] dy, 
                        double[:] dx2, double[:] dy2, double[:] dxdy):
    """
    Compute the Gaussian mode curvature values A,B,C at each (rx,ry) point.
    
    Inputs all have length N.
    
    Returns - (A, B, C, x', y', z')
                A,B,C are 1d arrays with the curvature coefficients. 
                x', y' and z' are (N,3) arrays giving the local basis vectors for the modes.
    """
    cdef:
        unsigned int N=len(rx), i
        
        np_.ndarray _A = np.empty((N), 'd')
        np_.ndarray _B = np.empty((N), 'd')
        np_.ndarray _C = np.empty((N), 'd')
        np_.ndarray _x = np.empty((N,3), 'd')
        np_.ndarray _y = np.empty((N,3), 'd')
        np_.ndarray _z = np.empty((N,3), 'd')
        
        double[:] A=_A,B=_B,C=_C
        double[:,:] x=_x,y=_y,z=_z
        
        vector_t a,b,c
        
    for i in range(N):
        #z' vector
        c.x = -dx[i]
        c.y = -dy[i]
        c.z = 1.0
        c = norm_(c)
        
        b.x = 0.0
        b.y = 0.0
        b.z = 1.0
        
        a = norm_(cross_(b,c))
        b = cross_(c,a)
        
        x[i,0] = a.x
        x[i,1] = a.y
        x[i,2] = a.z
        y[i,0] = b.x
        y[i,1] = b.y
        y[i,2] = b.z
        z[i,0] = c.x
        z[i,1] = c.y
        z[i,2] = c.z
        
        Fz = c.x*dx[i] + c.y*dy[i] + c.z
        
        ### This is d2F/dx^2
        A[i] = -( a.x*a.x*dx2[i] + 2*a.x*a.y*dxdy[i] + a.y*a.y*dy2[i] )/Fz
        ### This is d2F/dxdy
        B[i] = -( a.x*b.x*dx2[i] + (a.x*b.y + b.x*a.y)*dxdy[i] + a.y*b.y*dy2[i])/Fz
        ### This is d2F/dy^2
        C[i] = -( b.x*b.x*dx2[i] + 2*b.x*b.y*dxdy[i] + b.y*b.y*dy2[i] )/Fz
        
    return (_A, _B, _C, _x, _y, _z)
    
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
def build_interaction_matrix(double[:] rx, double[:] ry,
                             double[:] tx, double[:] ty,
                             double[:] A, double[:] B, double[:] C, 
                             double[:,:] x, double[:,:] y, double[:,:] z,
                             double wavelength, double spacing, 
                             double max_spacing, double blending,
                             double test_spacing):
    """
    Calculate the sparse matrix (COO format) to represent the complex amplitude of the field E_ij
    at mode i due to the mode j. E_ij = 1 where i==j.
    
    We calculate this in the local coordinate system of the decomposition plane.
    
    
    
    Parameters
    ----------
    
    rx - double[:] array of x-axis coordinates for the ray origins
    ry - double[:] array of y-axis coordinates for the ray origins
    A - double[:] wavefront curvature d2z'/dx'2 in ray-local basis
    B - double[:] wavefront curvature d2z'/dx'dy' in ray-local basis
    C - double[:] wavefront curvature d2z'/dy'2 in ray-local basis
    x,y,z - double[:,3] arrays giving the basis vectors for each ray (z is ray direction vector)
    wavelength - wavelength in microns
    spacing - the spacing between adjacent rays, in mm
    max_spacing - the maximum spacing between rays to calculate the cross-interaction. Controls the sparsity 
                    of the final result
    blending - adjusts the widths of each mode, relative to the spacing. width = spacing/blending
    test_spacing - the spacing used by the test-points
    
    Returns
    -------
    The returned tuple of arrays (data[:], (i[;],j[:]) is suitable to pass to scipy.sparse.coo_matrix().
    """
    
    cdef:
        unsigned int M, N=len(rx), L=len(tx), i, j, ct=0
        vector_t a,b,c,o,pt
        double complex[:] data
        unsigned long[:] xi, xj
        double phase=0.0, k, inv_root_area, impart
        double complex _A, _B, _C, detG0
        
    ##Estimate the size of the output array
    M = N*( 1 + int((7 * max_spacing * max_spacing)/(test_spacing*test_spacing)) ) 
    
    data = np.zeros(M, np.complex128)
    xi = np.zeros(M, np.uint)
    xj = np.zeros(M, np.uint)
    
    k = 2000.0*M_PI/wavelength
    
    impart = 2*(blending*blending)/(k*spacing*spacing)
    
    ###normalisation factor 1/root(area). Yes, there really is a double square root.
    inv_root_area = 1.0 #All modes have the same width
        
    for i in range(N):
        ### Origin of i'th mode
        o.x = rx[i]
        o.y = ry[i]
        o.z = 0.0
        
        a.x = x[i,0]
        a.y = x[i,1]
        a.z = x[i,2]
        b.x = y[i,0]
        b.y = y[i,1]
        b.z = y[i,2]
        c.x = z[i,0]
        c.y = z[i,1]
        c.z = z[i,2]
        
        ### Real parts are curvature, imaginary parts are 1/e^2 widths
        _A.real = A[i]
        _B.real = B[i]
        _C.real = C[i]
        
        _A.imag = impart
        _B.imag = 0 #It's a symmetric mode
        _C.imag = impart 
        
        detG0 = (_A*_C) #- (B*B)
        
        for j in range(L):
            pt.x = tx[j]
            pt.y = ty[j]
            pt.z = 0.0
            pt = subvv_(pt, o)
            
            if ((pt.x*pt.x) + (pt.y*pt.y)) <= (max_spacing*max_spacing):
                xi[ct] = i
                xj[ct] = j
                #if i == j:
                #    data[ct] = 1.0
                #else:
                data[ct] = calc_mode_U(_A, _B, _C, detG0, pt, a, b, c, k, phase, inv_root_area)
                ct += 1
    
    return (np.asarray(data)[:ct], (np.asarray(xi)[:ct], np.asarray(xj)[:ct]))


def apply_mode_curvature(GaussletCollection gc, double[:] A, double[:] B, double[:] C):
    """
    Adds the given wavefront curvature coefficients (i.e. 2nd derivatives) to the 
    given GaussletCollection object.
    """
    cdef:
        unsigned int i,j
        gausslet_t *g
        para_t *para
        ray_t *ray
        vector_t x,y,z, shift, r
        double px, py
        
    for i in range(gc.n_rays):
        g = &(gc.rays[i])
        ray = &(g.base_ray)
        z = ray.direction
        x = ray.E_vector
        y = cross_(x,z)
        for j in range(6):
            para = &(g.para[j])
            r = subvv_(para.origin, ray.origin)
            px = dotprod_(r, x)
            py = dotprod_(r, y)
            shift = multvs_(x, (-A[i]*px - B[i]*py))
            shift = addvv_(shift, multvs_(y, (-C[i]*py - B[i]*px)) )
            para.direction = norm_(addvv_(shift, para.direction))
            ###Not finished
    return 

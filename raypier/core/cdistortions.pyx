
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x ) nogil
    double pow(double x, double y) nogil
    double fabs(double) nogil
    double cos(double) nogil
    double sin(double) nogil
    double acos(double) nogil
    int isnan(double) nogil
    int abs(int) nogil

cdef extern from "float.h":
    double DBL_MAX
    
from libc.stdlib cimport malloc, free, realloc

cdef double INF=(DBL_MAX+DBL_MAX)
cdef double NAN=float("NaN")

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, Shape, Distortion
                
cimport cython
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
    
    
cdef packed struct zernike_coef_t:
    int j
    int n
    int m
    double value
    
    
cdef struct zer_n_m:
    int n
    int m
    

cdef zer_n_m eval_n_m(int j):
    ### Figure out n,m indices from a given j-index
    cdef:
        zer_n_m out
        
    return out

@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
cdef double eval_zernike_R(double r, int j, int n, int m, double[:] workspace) nogil:
    cdef:
        int jA, jB, jC, nA, nB, nC, mA, mB, mC
        double val
        
    if n < m:
        return 0.0
    
    val = workspace[j]
    if isnan(val):
        nA = n-1
        mA = abs(m-1)
        jA = (nA*(nA+2) + mA)//2
        nB = nA
        mB = m+1
        jB = (nB*(nB+2) + mB)//2
        val = r*(eval_zernike_R(r, jA, nA, mA, workspace) + eval_zernike_R(r, jB, nB, mB, workspace))
        nC = n-2
        mC = m
        jC = (nC*(nC+2) + mC)//2
        val -= eval_zernike_R(r, jC, nC, mC, workspace)
        workspace[j] = val
        return val
    else:
        return val
    
                
cdef class ZernikeDistortion(Distortion):
    cdef:
        public double unit_radius
        
        public int n_coefs, j_max
        
        zernike_coef_t[:] coefs
        double[:] workspace
         
    def __cinit__(self, **coefs):
        cdef:
            int n,m,i, j,size, j_max
            list clist
            zer_n_m nm
            double v
            
        clist = [(int(k[1:]), float(v)) for k,v in coefs if k.startswith("j")]
        j_max = max(a[0] for a in clist)
        self.j_max = j_max
        size = len(clist)
        self.coefs = <zernike_coef_t[:size]>malloc(size*sizeof(zernike_coef_t))
        self.n_coefs = size
        for i in range(size):
            j,v = clist[i]
            self.coefs[i].j = j
            nm = eval_n_m(j)
            self.coefs[i].n = nm.n
            self.coefs[i].m = nm.m
            self.coefs[i].value = v
            
        self.workspace = <double[:j_max]>malloc(j_max*sizeof(double))
            
    def __dealloc__(self):
        free(&(self.coefs[0]))
        free(&(self.workspace[0]))
        
    def __getitem__(self, int idx):
        if idx >= self.n_coefs:
            raise IndexError(f"Index {idx} greater than number of coeffs ({self.n_coefs}).")
        
    def __len__(self):
        return self.n_coefs
        
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    cdef double z_offset_c(self, double x, double y) nogil:
        cdef:
            double r, theta, Z=0.0, N, PH
            int i
            zernike_coef_t coef
            double[:] workspace
            
        x /= self.unit_radius
        y /= self.unit_radius
        r = sqrt(x*x + y*y)
        theta = atan2(y, x)
        
        workspace = self.workspace
        for i in range(self.j_max):
            workspace[i] = NAN
        workspace[0] = 1.0 #initialise the R_0_0 term
        
        for i in range(self.n_coefs):
            coef = self.coefs[i]
            
            if coef.m==0:
                N = sqrt(coef.n+1)
            else:
                N = sqrt(2*(coef.n+1))
                
            if coef.m >= 0:
                PH = cos(coef.m*theta)
            else:
                PH = -sin(coef.m*theta)
            
            Z += N * eval_zernike_R(r, coef.j, coef.n, abs(coef.m), workspace) * PH * coef.value
        
        return Z
        
        
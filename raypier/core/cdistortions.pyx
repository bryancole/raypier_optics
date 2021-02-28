
cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x ) nogil
    double pow(double x, double y) nogil
    double fabs(double) nogil
    double cos(double) nogil
    double sin(double) nogil
    double acos(double) nogil
    double floor(double) nogil
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
from cython.parallel cimport threadid
cimport numpy as np_
import numpy as np
cimport openmp


num_threads = openmp.omp_get_num_procs()
openmp.omp_set_num_threads(num_threads)
                
                
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
    int k
    double value
    
zcoefs_dtype = np.dtype([('j', np.int32),('n',np.int32),('m',np.int32),
                         ('k',np.int32),('value',np.double)])
    
    
cdef struct zer_nmk:
    int n
    int m
    int k
    
    
jnm_map = [(0,0,0)]
    

cdef zer_nmk eval_nmk_c(int j):
    ### Figure out n,m indices from a given j-index
    cdef:
        zer_nmk out
        int n, m, k, next_j, half_n
        
    if j<0:
        raise ValueError("J must be non-negative")
        
    if j >= len(jnm_map):
        n,m,k = jnm_map[-1]
        while True:
            while m<n:
                m += 2
                next_j = (n*(n+2) + m)//2
                half_n = int(floor(n/2.))
                k = half_n*(half_n + 1) + abs(m)
                jnm_map.append( (n,m,k) )
                if next_j == j:
                    out.n = n
                    out.m = m
                    out.k = k
                    return out
                elif next_j > j:
                    raise ValueError("Something went wrong here! Missed j value.")
            n += 1
            m = -2 - n
        
    nmk = jnm_map[j]
    out.n = nmk[0]
    out.m = nmk[1]
    out.k = nmk[2]
    return out


def eval_nmk(int j):
    cdef zer_nmk out
    out = eval_nmk_c(j)
    return (out.n, out.m, out.k)


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double zernike_R_c(double r, int k, int n, int m, double[:,:] workspace) nogil:
    cdef:
        int kA, kB, kC, nA, nB, nC, mA, mB, mC , half_n
        double val
        
    if n < m:
        return 0.0
    
    if n == 0:
        return 1.0
    
    val = workspace[0,k]
    if isnan(val):
        nA = n-1
        mA = abs(m-1)
        half_n = nA/2
        kA = half_n*(half_n+1) + abs(mA)#(nA*(nA+2) + mA)//2
        nB = nA
        mB = m+1
        kB = half_n*(half_n+1) + abs(mB) #(nB*(nB+2) + mB)//2
        nC = n-2
        mC = m
        half_n = nC/2
        kC = half_n*(half_n+1) + abs(mC) #(nC*(nC+2) + mC)//2
        val = r*(zernike_R_c(r, kA, nA, mA, workspace) + zernike_R_c(r, kB, nB, mB, workspace))
        val -= zernike_R_c(r, kC, nC, mC, workspace)
        workspace[0,k] = val
        return val
    else:
        return val

    
def zernike_R(double r, int k, int n, int m, double[:,:] workspace):
    return zernike_R_c(r,k,n,m, workspace)
    
    
@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double zernike_Rprime_c(double r, int k, int n, int m, double[:,:] workspace) nogil:
    cdef:
        int kA, kB, kC, nA, nB, nC, mA, mB, mC, half_n
        double val
        
    if n < m:
        return 0.0
    
    if n == 0:
        return 0.0
    
    val = workspace[1,k]
    if isnan(val):
        nA = n-1
        mA = abs(m-1)
        half_n = nA/2
        kA = half_n*(half_n+1) + abs(mA)
        nB = nA
        mB = m+1
        kB = half_n*(half_n+1) + abs(mB)
        nC = n-2
        mC = m
        kC = (half_n-1)*half_n + abs(mC)
        val = zernike_R_c(r, kA, nA, mA, workspace) + zernike_R_c(r, kB, nB, mB, workspace)
        val += r*( zernike_Rprime_c(r, kA, nA, mA, workspace) + \
                   zernike_Rprime_c(r, kB, nB, mB, workspace) )
        val -= zernike_Rprime_c(r, kC, nC, mC, workspace)
        workspace[1,k] = val
        return val
    else:
        return val
    
    
def zernike_Rprime(double r, int k, int n, int m, double[:,:] workspace):
    return zernike_Rprime_c(r,k,n,m,workspace)


# @cython.boundscheck(False)  # Deactivate bounds checking
# @cython.wraparound(False)   # Deactivate negative indexing.
# @cython.nonecheck(False)
# @cython.cdivision(True)
# cdef double zernike_R_over_r_c(double r, int k, int n, int m, double[:] workspace) nogil:
#     cdef:
#         int kA, kB, kC, nA, nB, nC, mA, mB, mC , half_n
#         double val, val1, val2, t1,t2,t3,t4
#         
#     if n < m:
#         return 0.0
#     
#     if (n == 0) and (r != 0.0):
#         return 1./r
# 
#     nA = n-1
#     mA = abs(m-1)
#     half_n = nA/2
#     kA = half_n*(half_n+1) + abs(mA)#(nA*(nA+2) + mA)//2
#     nB = nA
#     mB = m+1
#     kB = half_n*(half_n+1) + abs(mB) #(nB*(nB+2) + mB)//2
#     nC = n-2
#     mC = m
#     kC = (half_n-1)*half_n + abs(mC) #(nC*(nC+2) + mC)//2
#     
#     t1 = zernike_R_c(r, kA, nA, mA, workspace)
#     t2 = zernike_R_c(r, kB, nB, mB, workspace)
#     val1 = (t1 + t2)
#     
#     t3 = 0.0        
#     if (nA == 0) and (mA == 0):
#         t3 += 1.0
#     elif (nA == 2) and (mA == 0):
#         t3 -= 1.0
#     elif nA>0:
#         t3 += r*zernike_R_over_r_c(r, kA, nA, mA, workspace)
#     val2 = t3
#             
#     t4 = 0.0
#     if (nB == 0) and (mB == 0):
#         t4 += 1.0
#     elif (nB == 2) and (mB == 0):
#         t4 -= 1.0
#     elif nB>0:
#         t4 += r*zernike_R_over_r_c(r, kB, nB, mB, workspace)
#     val2 += t4
#     
#     val = val1
#     
#     val -= zernike_R_over_r_c(r, kC, nC, mC, workspace)
#         
#     with gil:
#         print("nmk:", (n,m,k), (nA,mA,kA),(nB,mB,kB), val1, val2, t1, t2, t3, t4 )
#     
#     return val


@cython.boundscheck(False)  # Deactivate bounds checking
@cython.wraparound(False)   # Deactivate negative indexing.
@cython.nonecheck(False)
@cython.cdivision(True)
cdef double zernike_R_over_r_c(double r, int k, int n, int m, double[:,:] workspace) nogil:
    """
    Results is undefined for m==0.
    """
    cdef:
        int kA, kB, kC, nA, nB, nC, mA, mB, mC , half_n
        double val
        
    if n < m:
        return 0.0
    
    val = workspace[2,k]
    if isnan(val):
        nA = n-1
        mA = abs(m-1)
        half_n = nA/2
        kA = half_n*(half_n+1) + abs(mA)
        nB = nA
        mB = m+1
        kB = half_n*(half_n+1) + abs(mB)
        nC = n-2
        mC = m
        half_n = nC/2
        kC = half_n*(half_n+1) + abs(mC)
        
        val = zernike_R_c(r, kA, nA, mA, workspace) + zernike_R_c(r, kB, nB, mB, workspace) 
        val -= zernike_R_over_r_c(r, kC, nC, mC, workspace)
        workspace[2,k]=val
    return val


def zernike_R_over_r(double r, int k, int n, int m, double[:,:] workspace):
    return zernike_R_over_r_c(r,k,n,m, workspace)

                
cdef class ZernikeDistortion(Distortion):
    cdef:
        public double unit_radius
        
        public int n_coefs, j_max, k_max
        
        public dict coef_map
        zernike_coef_t[:] coefs
        double[:,:,:] workspace
         
    def __cinit__(self, *args, **coefs):
        cdef:
            int n,m,i, j,size, j_max, k_max, half_n, n_max=0
            list clist
            zer_nmk nmk
            double v
            
        self.coef_map = {}
        self.unit_radius = coefs.get("unit_radius", 1.0)
        cdict = {int(k[1:]): float(v) for k,v in coefs.items() if k.startswith("j")}
        if args:
            for k,v in args[0]:
                cdict[int(k)] = float(v)
        self.k_max = 0
        self.set_coefs(list(cdict.items()))
        
    def __getitem__(self, int idx):
        cdef zernike_coef_t out
            
        if idx >= self.n_coefs:
            raise IndexError(f"Index {idx} greater than number of coeffs ({self.n_coefs}).")
        
        out = self.coefs[idx]
        return (out.j, out.n, out.m, out.k, out.value)
        
    def __len__(self):
        return self.n_coefs
    
    property workspace:
        def __get__(self):
            return np.asarray(self.workspace)
        
    def set_coefs(self, coefs):
        cdef:
            int i,j, n_max, k_max
            zernike_coef_t coef
            zer_nmk nmk
            double v
            
        clist = sorted(coefs)
        j_max = max(j for j,v in clist)
        self.j_max = j_max
        self.n_coefs = len(clist)
        self.coefs = np.empty(len(clist), dtype=zcoefs_dtype)
        n_max = 0
        self.coef_map = {}
        for i in range(self.n_coefs):
            j,v = clist[i]
            nmk = eval_nmk_c(j)
            self.coef_map[j]=i
            coef.j = j
            coef.n = nmk.n
            coef.m = nmk.m
            coef.k = nmk.k
            coef.value = v
            self.coefs[i] = coef
            if nmk.n > n_max:
                n_max = nmk.n
        k_max = (n_max//2)*(n_max//2 + 1) + n_max + 1
        if k_max > self.k_max:
            self.workspace = np.empty((num_threads, 3, k_max),'d')
            self.k_max = k_max
    
    def update_coef(self, int j, double value):
        cdef: 
            int i
        try:
            i = self.coef_map[j]
            self.coefs[i].value = value
        except KeyError:
            self.append_coef(j,value)
        
    def get_coef(self, int j):
        cdef int i=self.coef_map[j]
        return self.coefs[i].value
    
    def append_coef(self, int k, double value):
        raise NotImplementedError()
        
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cdef double z_offset_c(self, double x, double y) nogil:
        cdef:
            double r, theta, Z=0.0, N, PH, R
            int i
            zernike_coef_t coef
            double[:,:] workspace
            
        x /= self.unit_radius
        y /= self.unit_radius
        r = sqrt(x*x + y*y)
        theta = atan2(y, x)
        
        workspace = self.workspace[threadid(),:,:]
        for i in range(self.k_max):
            workspace[0,i] = NAN
        workspace[0,0] = 1.0 #initialise the R_0_0 term
        
        for i in range(self.n_coefs):
            coef = self.coefs[i]
            
            if coef.m==0:
                N = sqrt(coef.n+1)
            else:
                N = sqrt(2*(coef.n+1))
                
            N *= coef.value
                
            if coef.m >= 0:
                PH = cos(coef.m*theta)
            else:
                PH = -sin(coef.m*theta)
            
            R = zernike_R_c(r, coef.k, coef.n, abs(coef.m), workspace)
            Z += N * R * PH 
        return Z
        
    @cython.boundscheck(False)  # Deactivate bounds checking
    @cython.wraparound(False)   # Deactivate negative indexing.
    @cython.nonecheck(False)
    @cython.cdivision(True)
    cdef vector_t z_offset_and_gradient_c(self, double x, double y) nogil:
        """The z-axis surface sag is returned as the z-component 
        of the output vector. The x- and y-components of the surface
        gradient are placed in the x- and y- components of vector.
        I.e. the returned vector is
            (dz(x,y)/dx, dz(x,y)/dy, z(x,y))
        """
        cdef:
            double r, theta, N, PH, R, Rprime, R_over_r
            int i
            int thd = threadid()
            zernike_coef_t coef
            double[:,:] workspace
            vector_t Z
            
        x /= self.unit_radius
        y /= self.unit_radius
        r = sqrt(x*x + y*y)
        theta = atan2(y, x)
        
        workspace = self.workspace[thd,:,:]
        for i in range(self.k_max):
            workspace[0,i] = NAN
            workspace[1,i] = NAN
            workspace[2,i] = NAN
        
        Z.x = 0.0
        Z.y = 0.0
        Z.z = 0.0
        
        for i in range(self.n_coefs):
            coef = self.coefs[i]
                
            if coef.m >= 0:
                PH = cos(coef.m*theta)
                PHprime = -coef.m*sin(coef.m*theta)
            else:
                PH = -sin(coef.m*theta)
                PHprime = -coef.m*cos(coef.m*theta)
            
            R = zernike_R_c(r, coef.k, coef.n, abs(coef.m), workspace)
            Rprime = zernike_Rprime_c(r, coef.k, coef.n, abs(coef.m), workspace)
            R_over_r = zernike_R_over_r_c(r, coef.k, coef.n, abs(coef.m), workspace)
            
            if coef.m==0:
                N = sqrt(coef.n+1)
            else:
                N = sqrt(2*(coef.n+1))
            N *= coef.value
            
            Z.z += N * R * PH
            ### Now eval partial derivatives
            Z.x += N*(Rprime*cos(theta)*PH + R_over_r*(-sin(theta))*PHprime)
            Z.y += N*(Rprime*sin(theta)*PH + R_over_r*(cos(theta))*PHprime)
            
        Z.x /= self.unit_radius
        Z.y /= self.unit_radius
        
        return Z
        
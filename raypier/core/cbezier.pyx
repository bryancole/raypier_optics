
import numpy as np
cimport numpy as cnp

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c

from numpy import math


factorial = math.factorial


cdef double clip(double a, double lower, double upper) nogil:
    if lower < upper:
        if a < lower:
            a = lower
        if a > upper:
            a = upper
    else:
        if a > lower:
            a = lower
        if a < upper:
            a = upper
    return a


cdef double[:] binomial(int n):
    cdef int i
    return np.array([factorial(n)/(factorial(i)*factorial(n-i)) for i in range(n+1)])


cdef class BezierPatch():
    cdef:
        readonly int order_n, order_m
        double[:,:,:] _control_pts
        double[:] binom_n
        double[:] binom_m
        
        
    def __init__(self, unsigned int N, unsigned int M):
        self.order_n = <int>N
        self.order_m = <int>M
        self._control_pts = np.zeros((N+1,M+1,3))
        self.binom_n = binomial(N)
        self.binom_m = binomial(M)
        
    property control_pts:
        def __get__(self):
            return np.asarray(self._control_pts)
        
        def __set__(self, double[:,:,:] ctrl_pts):
            if ctrl_pts.shape != (self.order_n, self.order_m, 3):
                raise ValueError("Control points array must have shape (%d,%d,3)."%(self.order_n, self.order_m))
            self._control_pts = ctrl_pts.copy()
        
    cdef vector_t _eval_pt(self, double u, double v):
        cdef:
            int i,j, N=self.order_n, M=self.order_m
            vector_t out
            double[:] B_N = self.binom_n
            double[:] B_M = self.binom_m
            double coef
            double[:,:,:] ctrl=self._control_pts
            
        out.x = 0.0
        out.y = 0.0
        out.z = 0.0
            
        for i in range(self.order_n+1):
            for j in range(self.order_m+1):
                coef = B_N[i]*(u**i)*((1-u)**(N-i)) * B_M[j]*(v**j)*((1-v)**(M-j))
                out.x += coef*ctrl[i,j,0]
                out.y += coef*ctrl[i,j,1]
                out.z += coef*ctrl[i,j,2]
                
        return out
    
    def eval_pt(self, double u, double v):
        cdef:
            vector_t out
        out = self._eval_pt(u,v)
        return (out.x, out.y, out.z)
    
    def get_mesh(self, unsigned int N, unsigned int M):
        """
        Evaluates a triangular mesh over the patch area with resolution (N,M)
        for the two u,v coordinates
        
        returns : tuple(double[:,:] points, int[:,:] triangles)
        """
        
    
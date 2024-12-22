
#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

import numpy as np
cimport numpy as cnp

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c

from numpy import math


factorial = math.factorial


cdef struct vector3_t:
    vector_t p, dpdu, dpdv


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
        
        def __set__(self, pts_in):
            cdef:
                double[:,:,:] ctrl_pts=pts_in
            print("Got shape", ctrl_pts.shape, pts_in.shape)
            if pts_in.shape != (self.order_n+1, self.order_m+1, 3):
                raise ValueError("Control points array must have shape (%d,%d,3). Got %r"%(self.order_n+1, self.order_m+1, ctrl_pts.shape))
            self._control_pts = ctrl_pts #.copy()
        
    cdef vector_t _eval_pt(self, double u, double v) nogil:
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
    
    cdef vector3_t _eval_pt_and_grads(self, double u, double v) nogil:
        cdef:
            int i,j, N=self.order_n, M=self.order_m
            vector3_t out = [[0,0,0],[0,0,0],[0,0,0]]
            double[:] B_N = self.binom_n
            double[:] B_M = self.binom_m
            double coef, u_term, v_term
            double[:,:,:] ctrl=self._control_pts
            
            
        for i in range(self.order_n+1):
            for j in range(self.order_m+1):
                u_term = B_N[i]*(u**i)*((1-u)**(N-i))
                v_term = B_M[j]*(v**j)*((1-v)**(M-j))
                coef = u_term * v_term 
                out.p.x += coef*ctrl[i,j,0]
                out.p.y += coef*ctrl[i,j,1]
                out.p.z += coef*ctrl[i,j,2]
                
                #Take differential w.r.t. u
                coef = B_N[i]*(i-N*u)*(u**(i-1))*((1-u)**(N-1-i)) * v_term
                out.dpdu.x += coef*ctrl[i,j,0]
                out.dpdu.y += coef*ctrl[i,j,1]
                out.dpdu.z += coef*ctrl[i,j,2]
                
                #Take differential w.r.t. v
                coef = u_term * B_M[j]*(j-M*v)*(v**(j-1))*((1-v)**(M-1-j))
                out.dpdv.x += coef*ctrl[i,j,0]
                out.dpdv.y += coef*ctrl[i,j,1]
                out.dpdv.z += coef*ctrl[i,j,2]
                
        return out
    
    def eval_pt(self, double u, double v):
        cdef:
            vector_t out
        out = self._eval_pt(u,v)
        return (out.x, out.y, out.z)
    
    def eval_pt_and_grad(self, double u, double v):
        cdef:
            int i
            vector3_t out
        out = self._eval_pt_and_grad(u,v)
        return [(out.p.x, out.p.y, out.p.z),
                (out.dpdu.x, out.dpdu.y, out.dpdu.z),
                (out.dpdv.x, out.dpdv.y, out.dpdv.z)]
    
    def get_mesh(self, unsigned int N, unsigned int M):
        """
        Evaluates a triangular mesh over the patch area with resolution (N,M)
        for the two u,v coordinates
        
        returns : tuple(double[:,3] points, int[:,3] triangles, double[:,2] uv_coords)
        """
        cdef:
            double du=1./(N-1), dv=1./(M-1), u, v
            vector_t pt
            double[:,:,:] out
            double[:,:,:] uv_out
            long[:,:] pt_ids
            long long[:,:] cells
            int i,j, ct=0
            
        points = np.empty((N,M,3), dtype=np.float64)
        uv = np.empty((N,M,2), dtype=np.float64)
        out = points
        uv_out = uv
        pt_ids = np.arange(N*M).reshape(N,M)
        cells = np.empty(((N-1)*(M-1)*2,3), np.int64)
        
        with nogil:
            for i in range(N):
                for j in range(M):
                    u = i*du
                    v = j*dv
                    uv_out[i,j,0] = u
                    uv_out[i,j,1] = v
                    pt = self._eval_pt(u, v)
                    out[i,j,0] = pt.x
                    out[i,j,1] = pt.y
                    out[i,j,2] = pt.z
                    
            for i in range(N-1):
                for j in range(M-1):
                    cells[ct,0] = pt_ids[i,j]
                    cells[ct,1] = pt_ids[i,j+1]
                    cells[ct,2] = pt_ids[i+1,j+1]
                    ct += 1
                    cells[ct,0] = pt_ids[i,j]
                    cells[ct,1] = pt_ids[i+1,j+1]
                    cells[ct,2] = pt_ids[i+1,j]
                    ct += 1
                    
        return (points.reshape(-1,3), np.asarray(cells), uv.reshape(-1,2))
    
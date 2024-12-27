
#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

import numpy as np
cimport numpy as cnp

from .ctracer cimport Face, sep_, \
        vector_t, ray_t, FaceList, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, mag_, transform_t, Transform, transform_c,\
                rotate_c, intersect_t, cross_, uv_coord_t, invert_
                
from .obbtree cimport OBBTree, intersection

from numpy import math


cdef intersect_t NO_INTERSECTION

NO_INTERSECTION.dist = -1
NO_INTERSECTION.face_idx = -1
NO_INTERSECTION.piece_idx = 0


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
    
    
cdef class BezierPatchFace(Face):
    cdef:
        public BezierPatch bezier
        public OBBTree obbtree
        long long[:] workspace
        #double[:,:] normals #I'm not sure I need normals
        
        ### Set the resolution of the mesh
        unsigned int u_res, v_res
        double[:,:] uvs
        public double atol
        public int invert_normals
        
    def __cinit__(self, **kwds):
        self.u_res = kwds.get("u_res", 20)
        self.v_res = kwds.get("v_res", 20)
        self.atol = kwds.get("atol", 1e-10)
        self.invert_normals = 1 if kwds.get("invert_normals", 0) else 0
        bezier = kwds["bezier"]
        self.bezier = bezier
        pts, cells, uvs = bezier.get_mesh(self.u_res, self.v_res)
        
        tree = OBBTree(pts.copy(), np.ascontiguousarray(cells, dtype=np.int32))
        tree.max_level = kwds.get("max_level",100)
        tree.number_of_cells_per_node = kwds.get("cells_per_node", 2)
        if tree.level <= 0:
            tree.build_tree()
        self.obbtree = tree
        self.workspace = np.zeros(tree.level+1, np.int64)
        
        self.uvs = uvs
        
    cdef uv_coord_t interpolate_cell(self, int cell_idx, vector_t pt):
        cdef:
            OBBTree tree=self.obbtree
            int[:] cell
            double x2,x3,y3, alpha0, alpha1, alpha2
            double[:,:] uvs=self.uvs
            vector_t p1, p2, p3, edge1, edge2, en1, en2
            uv_coord_t out
            
        ###Now find the UV-coord of the cell intersection by linear barycentric interpolation
        cell = tree.cells[cell_idx]
        p1 = tree.get_point_c(cell[0])
        p2 = tree.get_point_c(cell[1])
        p3 = tree.get_point_c(cell[2])
        
        pt = subvv_(pt,p1)
        
        edge1 = subvv_(p2,p1)
        edge2 = subvv_(p3,p1)
        
        en1 = norm_(edge1)
        en2 = norm_(cross_(en1, cross_(edge1, edge2)))
        
        x2 = mag_(edge1)
        x3 = dotprod_(edge2, en1)
        y3 = dotprod_(edge2, en2)
        
        #Barycentric coefficients
        alpha0 = -pt.x/x2 - pt.y/y3 + pt.y*x3/(x2*y3) + pt.z
        alpha1 = (pt.x*y3 - pt.y*x3)/(x2*y3)
        alpha2 = pt.y/y3
        print("alpha:", alpha0 + alpha1 + alpha2)
        out.u = alpha0 * uvs[cell[0],0] + alpha1 * uvs[cell[1],0] + alpha2 * uvs[cell[2],0]
        out.v = alpha0 * uvs[cell[0],1] + alpha1 * uvs[cell[1],1] + alpha2 * uvs[cell[2],1]
        return out
        
        
    cdef intersect_t intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
        """Intersects the given ray with this face.
        
        params:
          p1 - the origin of the input ray, in local coords
          p2 - the end-point of the input ray, in local coords
          
        returns:
          the distance along the ray to the first valid intersection. No
          intersection can be indicated by a negative value.
        """
        cdef:
            intersection it #OBBTree specific intersection
            intersect_t out #Face intersection data
            # cdef struct intersect_t:
            #     double dist
            #     int face_idx
            #     int piece_idx
            
            OBBTree tree=self.obbtree
            BezierPatch bezier=self.bezier
            int[:] cell
            double du, dv, dist, tol=self.atol*self.atol
            double[:,:] uvs=self.uvs
            vector_t pt, d, normal
            uv_coord_t uv
            vector3_t pt_grads
            
        it = tree.intersect_with_line_c(p1, p2, self.workspace)
        if it.alpha < self.tolerance:
            print("no mesh intersection")
            return NO_INTERSECTION
        
        d = norm_(subvv_(p2,p1))
        pt = addvv_(multvs_(p2, it.alpha), multvs_(p1, 1.0-it.alpha)) 
        
        #out.dist = it.alpha * mag_(subvv_(p2,p1))
        #out.piece_idx = <int>(it.cell_idx) #casting long to int. Hopefully there are not too many cells
        
        ###Now find the UV-coord of the cell intersection by linear barycentric interpolation
        # This doesn't work right!
        uv = self.interpolate_cell(it.cell_idx, pt)
        print("mesh intersection:", it.alpha, uv.u, uv.v)
        
        for i in range(100):
            pt_grads = bezier._eval_pt_and_grads(uv.u, uv.v)
            
            normal = norm_(cross_(pt_grads.dpdu, pt_grads.dpdv))
            dist = dotprod_(subvv_(pt_grads.p, p1), normal)/dotprod_(d,normal)
            dp = subvv_(addvv_(p1, multvs_(d, dist)), pt_grads.p)
            du = dotprod_(pt_grads.dpdu, dp) / mag_sq_(pt_grads.dpdu)
            dv = dotprod_(pt_grads.dpdv, dp) / mag_sq_(pt_grads.dpdv)
            
            uv.u += du
            uv.v += dv
            
            if (du*du < tol) and (dv*dv < tol):
                break
        else:
            return NO_INTERSECTION
        print("iteractions:", i, du, dv)
        pt = bezier._eval_pt(uv.u, uv.v)
        
        out.dist = mag_(subvv_(pt,p1))
        out.uv = uv
        return out
    
    cdef void compute_normal_and_tangent_c(self, vector_t p, intersect_t *it, vector_t *normal, vector_t *tangent):
        cdef:
            BezierPatch bezier=self.bezier
            uv_coord_t uv
            vector3_t pt_grads
            vector_t n
            
        uv = it.uv
        pt_grads = bezier._eval_pt_and_grads(uv.u, uv.v)
        
        n = norm_(cross_(pt_grads.dpdu, pt_grads.dpdv))
        normal[0] = invert_(n) if self.invert_normals else n
        tangent[0] = norm_(pt_grads.dpdu)
        
        
        

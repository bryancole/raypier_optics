
"""
Implementation of a Oriented Boundary Boxx Tree (OBB Tree) spatial search algorithm.
Somewhat copied from the vtkOBBTree implementation.
"""
from examples.entry_ap_demo import ratio
from test.test_obbtree import TestComputeOBB

cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    
import numpy as np
cimport numpy as np_
from libc.float cimport DBL_MAX, DBL_MIN
from cython cimport view

from .ctracer cimport vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_


cdef struct MaxElem:
    double Amax
    int i
    int j 
    

cdef vector_t as_vector(double[3] pt):
    cdef:
        vector_t v
    v.x = pt[0]
    v.y = pt[1]
    v.z = pt[2]
    return v
    

cdef MaxElem max_elem(double[:,:] A):
    """Find max off-diagonal element of square matrix A
    """
    cdef:
        int i,j, n=A.shape[0]
        double absA
        MaxElem out
        
    out.Amax = 0.0
    for i in range(n-1):
        for j in range(i+1,n):
            absA = fabs(A[i,j])
            if absA >= out.Amax: 
                out.Amax = absA
                out.i = i
                out.j = j
    return out

cdef rotate(double[:,:] A, double[:,:] p, int k, int l):
    """
    Rotate matrix A by transform p, such that A[k,l]==0
    """
    cdef:
        int n, i
        double Adiff, phi, c, tau, temp
        
    n = len(A)
    Adiff = A[l,l] - A[k,k]
    if fabs(A[k,l]) < fabs(Adiff)*1.0e-36: 
        t = A[k,l]/Adiff
    else:
        phi = Adiff/(2.0*A[k,l])
        t = 1.0/(abs(phi) + sqrt(phi**2 + 1.0))
        if phi < 0.0: 
            t = -t
    c = 1.0/sqrt(t**2 + 1.0)
    s = t*c
    tau = s/(1.0 + c)
    temp = A[k,l]
    A[k,l] = 0.0
    A[k,k] = A[k,k] - t*temp
    A[l,l] = A[l,l] + t*temp
    for i in range(k):
    # Case of i < k
        temp = A[i,k]
        A[i,k] = temp - s*(A[i,l] + tau*temp)
        A[i,l] = A[i,l] + s*(temp - tau*A[i,l])
    for i in range(k+1,l): # Case of k < i < l
        temp = A[k,i]
        A[k,i] = temp - s*(A[i,l] + tau*A[k,i])
        A[i,l] = A[i,l] + s*(temp - tau*A[i,l])
    for i in range(l+1,n): # Case of i > l
        temp = A[k,i]
        A[k,i] = temp - s*(A[l,i] + tau*temp)
        A[l,i] = A[l,i] + s*(temp - tau*A[l,i])
    for i in range(n):
    # Update transformation matrix
        temp = p[i,k]
        p[i,k] = temp - s*(p[i,l] + tau*p[i,k])
        p[i,l] = p[i,l] + s*(temp - tau*p[i,l])
        

cdef int _jacobi(double[:,:] A, double[:] eigenvals, double[:,:] p, double tol=1e-9):
    cdef:
        int i, j, n = A.shape[0]
        int max_rot = 5*(n**2)
        MaxElem me
        
    for i in range(n):
        for j in range(n):
            if i==j:
                p[i,j] = 1.0
            else:
                p[i,j] = 0.0
        
    for i in range(max_rot):
        me = max_elem(A)
        if me.Amax < tol:
            for j in range(n):
                eigenvals[j] = A[j,j]
            return 0
        rotate(A,p,me.i, me.j)
    return -1


cdef inline vector_t as_vect(double[:] v_in):
    cdef:
        vector_t out
    out.x = v_in[0]
    out.y = v_in[1]
    out.z = v_in[2]
    return out
    
    
    
def jacobi(double[:,:] A, double tol=1e-9):
    cdef:
        int n=len(A)
        double[:,:] p = np.identity(n)*1.0
        double[:] vals = np.empty(n, dtype='d')
        
    if _jacobi(A, vals, p, tol=tol)<0:
        raise Exception("Jacobi method did not converge")
    
    return (vals, np.asarray(p))

def jacobi2(double[:,:] A, double tol=1e-9):
    cdef:
        int i, n = len(A)
        int max_rot = 5*(n**2)
        double[:,:] p = np.identity(n)*1.0
        MaxElem me
        
    for i in range(max_rot):
        me = max_elem(A)
        if me.Amax < tol:
            return np.diagonal(A),np.asarray(p)
        rotate(A,p,me.i, me.j)
    raise Exception("Jacobi method did not converge")



cdef class OBBNode(object):
    cdef:
        vector_t _corner #one corner of the box
        double[3][3] _axes #the three primary axes of the box, ordered from longest to smallest
        public OBBNode parent, child1, child2 #children will be None for leaf nodes. Parent will be None for root
        public int[:] cell_list
        
    property corner:
        def __get__(self):
            return np.array([self._corner.x, self._corner.y, self._corner.z])
        
        def __set__(self, c):
            self._corner.x = c[0]
            self._corner.y = c[1]
            self._corner.z = c[2]
            
    property axes:
        def __get__(self):
            cdef:
                int i,j
                double[:,:] axes = np.empty((3,3), 'd')
            for i in range(3):
                for j in range(3):
                    axes[i,j] = self._axes[i][j]
            return np.asarray(axes)
        
        def __set__(self, double[:,:] axes):
            cdef:
                int i,j
            for i in range(3):
                for j in range(3):
                    self._axes[i][j] = axes[i,j]
        
        
cdef class OBBTree(object):
    cdef:
        double[:,:] points
        int[:,:] cells #always triangles
        int[:] point_mask
        OBBNode root
        
        public int max_level
        public int number_of_cells_per_node
        
    def __init__(self, double[:,:] points, int[:,:] cells):
        self.points = points
        self.cells = cells
        self.point_mask = np.zeros(points.shape[0], dtype=np.int32)
        self.max_level = 12
        self.number_of_cells_per_node = 1
        
    def build_tree(self):
        cdef:
            int n_cells = self.cells.shape[0]
            int[:] cell_ids = np.arange(n_cells)
            OBBNode obb

        obb = self.compute_obb_cells(cell_ids)
        self.root = self.c_build_tree(obb, 0)
        
    cdef c_build_tree(self, OBBNode obb, int level):
        cdef:
            OBBNode child1, child2
            int i, j, split_plane, icell
            vector_t p, n, c, pt
            double best_ratio, ratio
            int[:] cell
            int[:] cell_ids = obb.cell_list
            int split_acceptable=0, n_cells = cell_ids.shape[0]
            int negative, positive, out_cells_left=0, out_cells_right=(n_cells-1)
            view.array out_cells = view.array(shape=cell_ids.shape, 
                                              itemsize=sizeof(int),
                                              format="i")
        
        if (level < self.max_level) and (n_cells > self.number_of_cells_per_node):
            
            p.x = obb._corner.x + obb._axes[0][0]/2. + obb._axes[1][0]/2. + obb._axes[2][0]/2.
            p.y = obb._corner.y + obb._axes[0][1]/2. + obb._axes[1][1]/2. + obb._axes[2][1]/2.
            p.z = obb._corner.z + obb._axes[0][2]/2. + obb._axes[1][2]/2. + obb._axes[2][2]/2.
                 
            best_ratio = 1.0
            
            for split_plane in range(3):
                n = norm_(as_vector(obb._axes[split_plane]))
                out_cell_left=0
                out_cell_right=(n_cells-1)
                for i in range(n_cells):
                    cell = self.cells[i]
                    c.x = 0
                    c.y = 0
                    c.z = 0
                    negative = 0
                    positive = 0
                    for j in range(3): #Always triangles
                        icell = cell[j]
                        pt.x = self.points[icell, 0]
                        pt.y = self.points[icell, 1]
                        pt.z = self.points[icell, 2]
                        
                        val = dotprod_(n, subvv_(pt, p))
                        
                        c.x += pt.x
                        c.y += pt.y
                        c.z += pt.z
                        
                        if val < 0:
                            negative = 1
                        else:
                            positive = 1
                            
                    if (negative and positive):
                        c.x /= 3
                        c.y /= 3
                        c.z /= 3
                        if dotprod_(n, subvv_(c, p)) < 0:
                            out_cells[out_cells_left] = i
                            out_cells_left += 1
                        else:
                            out_cells[out_cells_right] = i
                            out_cells_right -= 1
                    else:
                        if negative:
                            out_cells[out_cells_left] = i
                            out_cells_left += 1
                        else:
                            out_cells[out_cells_right] = i
                            out_cells_right -= 1
                            
                ratio = fabs( <double>(n_cells - out_cells_right - out_cells_left - 1)/n_cells )
                
                if (ratio < 0.6):
                    split_acceptable = 1
                    break
                else:
                    if ratio < best_ratio:
                        best_ratio = ratio
                        best_plane = split_plane
                        
            if best_ratio < 0.95:
                split_plane = best_plane
                split_acceptable = 1
                
            if split_acceptable:
                child1 = self.compute_obb_cells(out_cells[:out_cells_left])
                child2 = self.compute_obb_cells(out_cells[out_cells_right+1:])
                obb.child1 = child1
                obb.child2 = child2
                child1.parent = obb
                child2.parent = obb
                self.c_build_tree(child1, level+1)
                self.c_build_tree(child2, level+1)
        
        
    def compute_obb(self, int[:] point_ids):
        """
        Compute an OBB from the list of points given. Return the corner point
        and the three axes defining the orientation of the OBB. Also return
        a sorted list of relative "sizes" of axes for comparison purposes.
        
        In fact, we want to compute the OBB based on cell 'moments'. These are 
        weighed by cell area.
        """
        cdef:
            double[3] pt0
            double[:] pt
            int i,j, n_pts = point_ids.size
            ###co-variance matrix
            double[3][3] a = [[0.0,0.0,0.0],
                              [0.0,0.0,0.0],
                              [0.0,0.0,0.0]]
            double[3] axis_lens
            OBBNode node = OBBNode()
            double[3] tmin = [DBL_MAX, DBL_MAX, DBL_MAX]
            double[3] tmax = [-DBL_MAX, -DBL_MAX, -DBL_MAX]
            double t
            vector_t ax, vpt, centre, corner
        
        centre.x=0.
        centre.y=0.
        centre.z=0.
        for i in range(n_pts):
            pt = self.points[point_ids[i]]
            centre.x += pt[0]
            centre.y += pt[1]
            centre.z += pt[2]
        centre.x /= n_pts
        centre.y /= n_pts
        centre.z /= n_pts

        ### Populate covariance matrix
        for i in range(n_pts):
            pt = self.points[point_ids[i]]
            pt0[0] = pt[0] - centre.x
            pt0[1] = pt[1] - centre.y
            pt0[2] = pt[2] - centre.z
            for j in range(3):
                a[0][j] += pt0[0] * pt0[j]
                a[1][j] += pt0[1] * pt0[j]
                a[2][j] += pt0[2] * pt0[j]
        for j in range(3):
            a[0][j] /= n_pts
            a[1][j] /= n_pts
            a[2][j] /= n_pts
            
        ### Find eigenvectors
        if _jacobi(a, <double[:]>axis_lens, node._axes) < 0:
            raise Exception("Jacobi iteration failed.")
        ### Jacobi method outputs unit-magnitude eigenvectors so we don't need to re-normalise
        
        ### Get bounding box by projecting onto eigenvectors
        for i in range(n_pts):
            vpt = as_vect(self.points[point_ids[i]])
            vpt = subvv_(vpt, centre)
            for j in range(3): #iterate over axes
                ax = as_vect(node.axes[j])
                t = dotprod_(ax, vpt)
                if t < tmin[j]:
                    tmin[j] = t
                if t > tmax[j]:
                    tmax[j] = t
        
        print("min: ", tmin[0], tmin[1], tmin[2], "    max: ", tmax[0], tmax[1], tmax[2])
        corner = addvv_(centre, multvs_(as_vect(node._axes[0]), tmin[0]))
        corner = addvv_(corner, multvs_(as_vect(node._axes[1]), tmin[1]))
        corner = addvv_(corner, multvs_(as_vect(node._axes[2]), tmin[2]))
            
        for j in range(3):
            for i in range(3):
                node._axes[j][i] *= (tmax[j] - tmin[j])
                #pt0 = node.axes[i]
                #pt0[j] *= (tmax[j] - tmin[j])
        
        node._corner = corner
        return node
    
    cdef clear_point_mask(self):
        cdef:
            int i, n=self.point_mask.shape[0]
        for i in range(n):
            self.point_mask[i] = 0
    
    def compute_obb_cells(self, int[:] cell_list):
        """
        cell_list - a 1d array of indices into the cells member.
        
        returns:
            - a OBBNode instance
        """
        cdef:
            vector_t mean, p, q, r, dp0, dp1, c
            int n_cells = cell_list.shape[0]
            int i,j, n_pts=self.points.shape[0]
            double tot_mass=0.0, tri_mass
            int[:] tri
            double[3] axis_lens
            double[3] _mean
            double[3] tmin = [DBL_MAX, DBL_MAX, DBL_MAX]
            double[3] tmax = [-DBL_MAX, -DBL_MAX, -DBL_MAX]
            ###co-variance matrix
            double[3][3] a = [[0.0,0.0,0.0],
                              [0.0,0.0,0.0],
                              [0.0,0.0,0.0]]
            double[:] a0 = a[0]
            double[:] a1 = a[1]
            double[:] a2 = a[2]
            OBBNode node = OBBNode()
            
        mean.x=0.0
        mean.y=0.0
        mean.z=0.0
        
        self.clear_point_mask()
            
        for i in range(n_cells):
            tri = self.cells[i]
            self.point_mask[tri[0]] = 1
            self.point_mask[tri[1]] = 1
            self.point_mask[tri[2]] = 1
            p = as_vect(self.points[tri[0]])
            q = as_vect(self.points[tri[1]])
            r = as_vect(self.points[tri[2]])
            
            dp0 = subvv_(q,p)
            dp1 = subvv_(r,p)
            
            c.x = (p.x + q.x + r.x)/3.
            c.y = (p.y + q.y + r.y)/3.
            c.z = (p.z + q.z + r.z)/3.
            
            tri_mass = 0.5 * dotprod_(dp0, dp1)
            tot_mass += tri_mass
            
            mean.x += tri_mass * c.x
            mean.y += tri_mass * c.y
            mean.z += tri_mass * c.z
            
            ## on-diagonal terms
            a0[0] += tri_mass * (9 * c.x * c.x + p.x * p.x + q.x * q.x + r.x * r.x) / 12;
            a1[1] += tri_mass * (9 * c.y * c.y + p.y * p.y + q.y * q.y + r.y * r.y) / 12;
            a2[2] += tri_mass * (9 * c.z * c.z + p.z * p.z + q.z * q.z + r.z * r.z) / 12;
            
            ## off-diagonal terms
            a0[1] += tri_mass * (9 * c.x * c.y + p.x * p.y + q.x * q.y + r.x * r.y) / 12;
            a0[2] += tri_mass * (9 * c.x * c.z + p.x * p.z + q.x * q.z + r.x * r.z) / 12;
            a1[2] += tri_mass * (9 * c.y * c.z + p.y * p.z + q.y * q.z + r.y * r.z) / 12;
            
        mean.x /= tot_mass
        mean.y /= tot_mass
        mean.z /= tot_mass
        
        ## matrix is symmetric
        a1[0] = a0[1];
        a2[0] = a0[2];
        a2[1] = a1[2];
        
        _mean[0] = mean.x
        _mean[1] = mean.y
        _mean[2] = mean.z
        
        ## get covariance from moments
        for i in range(3):
            for j in range(3):
                a[i][j] = (a[i][j] / tot_mass) - (_mean[i] * _mean[j])
                
        ### Find eigenvectors
        if _jacobi(a, <double[:]>axis_lens, node._axes) < 0:
            raise Exception("Jacobi iteration failed.")
        ### Jacobi method outputs unit-magnitude eigenvectors so we don't need to re-normalise
        
        ### Get bounding box by projecting onto eigenvectors
        for i in range(n_pts):
            if self.point_mask[i]:
                vpt = as_vect(self.points[i])
                vpt = subvv_(vpt, mean)
                for j in range(3): #iterate over axes
                    ax = as_vect(node.axes[j])
                    t = dotprod_(ax, vpt)
                    if t < tmin[j]:
                        tmin[j] = t
                    if t > tmax[j]:
                        tmax[j] = t
        
        #print("min: ", tmin[0], tmin[1], tmin[2], "    max: ", tmax[0], tmax[1], tmax[2])
        corner = addvv_(mean, multvs_(as_vect(node._axes[0]), tmin[0]))
        corner = addvv_(corner, multvs_(as_vect(node._axes[1]), tmin[1]))
        corner = addvv_(corner, multvs_(as_vect(node._axes[2]), tmin[2]))
            
        for j in range(3):
            for i in range(3):
                node._axes[j][i] *= (tmax[j] - tmin[j])
                #pt0 = node.axes[i]
                #pt0[j] *= (tmax[j] - tmin[j])
        
        node._corner = corner
        return node
        
            
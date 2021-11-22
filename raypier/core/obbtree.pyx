#cython: language_level=3
"""
Implementation of a Oriented Boundary Boxx Tree (OBB Tree) spatial search algorithm.
Somewhat copied from the vtkOBBTree implementation.
"""

cdef extern from "math.h":
    double M_PI
    double sqrt(double) nogil
    double atan2 (double y, double x )
    double pow(double x, double y)
    double fabs(double)
    
from libc.stdlib cimport malloc, free, realloc
    
import numpy as np
cimport numpy as np_
from libc.float cimport DBL_MAX, DBL_MIN
from cython cimport view, boundscheck, cdivision, initializedcheck

from .ctracer cimport vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_


cdef struct MaxElem:
    double Amax
    int i
    int j 
    

cdef inline vector_t as_vector(double[3] pt):
    cdef:
        vector_t v
    v.x = pt[0]
    v.y = pt[1]
    v.z = pt[2]
    return v
    
    
@boundscheck(False)
@initializedcheck(False)
@cdivision(True)
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


@boundscheck(False)
@initializedcheck(False)
@cdivision(True)
cdef void rotate(double[:,:] A, double[:,:] p, int k, int l):
    """
    Rotate matrix A by transform p, such that A[k,l]==0
    """
    cdef:
        int n, i
        double Adiff, phi, c, tau, temp, t
        
    n = A.shape[0] #len(A)
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
        
        
@boundscheck(False)
@initializedcheck(False)
@cdivision(True)
cdef int _jacobi(double[:,:] A, double[:] eigenvals, double[:,:] p, double tol):
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


@boundscheck(False)
@initializedcheck(False)
@cdivision(True)
cdef int _jacobi_no_eigenvals(double[:,:] A, double[:,:] p, double tol):
    cdef:
        int i, j, n = A.shape[0]
        int max_rot = 5*(n*n)
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
            return 0
        rotate(A,p,me.i, me.j)
    return -1


@boundscheck(False)
@initializedcheck(False)
@cdivision(True)
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



cdef struct obbnode_t:
    vector_t _corner #one corner of the box
    vector_t[3] _axes #the three primary axes of the box, ordered from longest to smallest
    long parent, child1, child2 #children will be -1 for leaf nodes. Parent will be -1 for root
    int cell_list_idx
    int n_cells



cdef class OBBNode(object):
    cdef:
        vector_t _corner #one corner of the box
        vector_t[3] _axes #the three primary axes of the box, ordered from longest to smallest
        public OBBNode parent, child1, child2 #children will be None for leaf nodes. Parent will be None for root
        public int cell_list_idx
        public int n_cells
        
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
                int i
            out = np.empty(shape=(3,3))
            for i in range(3):
                out[i,0] = self._axes[i].x
                out[i,1] = self._axes[i].y
                out[i,2] = self._axes[i].z
            return out
        
        def __set__(self, double[:,:] axes):
            cdef:
                int i
            for i in range(3):
                self._axes[i].x = axes[i,0]
                self._axes[i].y = axes[i,1]
                self._axes[i].z = axes[i,2]
                    
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef int intersect_line_c(self, vector_t p1, vector_t p2):
        cdef:
            int i
            double rangeAmin, rangeAmax, rangePmin, rangePmax, dotP
            vector_t ax
            
        for i in range(3):
            ax = self._axes[i]
            
            rangeAmin = dotprod_(self._corner, ax)
            rangeAmax = rangeAmin + mag_sq_(ax)
            
            rangePmin = dotprod_(p1, ax)
            rangePmax = rangePmin
            dotP = dotprod_(p2, ax)
            if dotP < rangePmin:
                rangePmin = dotP
            else:
                rangePmax = dotP
            
            if (rangeAmax < rangePmin) or (rangePmax < rangeAmin):
                return 0
        return 1
    
    def intersect_line(self, p1, p2):
        cdef:
            vector_t _p1, _p2
        _p1.x = p1[0]
        _p1.y = p1[1]
        _p1.z = p1[2]
        _p2.x = p2[0]
        _p2.y = p2[1]
        _p2.z = p2[2]
        return bool(self.intersect_line_c(_p1, _p2))            
                    
    def clear_children(self):
        ch1 = self.child1
        ch2 = self.child2
        if ch1 is not None:
            ch1.clear_children()
            del ch1.parent
            del self.child1
        if ch2 is not None:
            ch2.clear_children()
            del ch2.parent
            del self.child2
            
        
cdef class OBBTree(object):
    cdef:
        double[:,:] points
        int[:,:] cells #always triangles
        int[:] cell_ids
        int[:] point_mask
        public OBBNode root
        
        public int max_level
        public int number_of_cells_per_node
        
        double[:,:] _covar #A workspace for obb calculation
        double[:,:] _axes
        
        unsigned long n_nodes
        unsigned long max_n_nodes
        obbnode_t *nodes
    
    def __dealloc__(self):
        free(self.nodes)
        
    def __cinit__(self, double[:,:] points, int[:,:] cells):
        cdef:
            int max_n_nodes
            
        self.points = points
        self.cells = cells
        self.cell_ids = np.arange(cells.shape[0])
        
        self.point_mask = np.zeros(points.shape[0], dtype=np.int32)
        self.max_level = 12
        self.number_of_cells_per_node = 1
        self._covar = np.zeros(shape=(3,3))
        self._axes = np.zeros(shape=(3,3))
        
        max_n_nodes = 2**(int(np.ceil(np.log2(cells.shape[0]))) + 1)
        self.max_n_nodes = max_n_nodes
        self.nodes = <obbnode_t*>malloc(max_n_nodes*sizeof(obbnode_t))
        self.n_nodes = 0
        
    cdef obbnode_t get_node_c(self, unsigned long i):
        return self.nodes[i]
    
    cdef void set_node_c(self, unsigned long i, obbnode_t node):
        self.nodes[i] = node
        
    cdef unsigned long add_node_c(self, obbnode_t r):
        cdef:
            n_nodes = self.n_nodes
        if self.n_nodes == self.max_n_nodes:
            if self.max_n_nodes == 0:
                self.max_n_nodes = 1
            else:
                self.max_n_nodes *= 2
            self.nodes = <obbnode_t*>realloc(self.nodes, self.max_n_nodes*sizeof(obbnode_t))
        self.nodes[n_nodes] = r
        self.n_nodes += 1
        return n_nodes
        
    def clear_tree(self):
        self.n_nodes = 0
        
    def build_tree(self):
        cdef:
            int n_cells = self.cells.shape[0]
            obbnode_t obb

        self.n_nodes = 0
        obb = self.compute_obb_cells_c(0, n_cells)
        self.root = self.c_build_tree(obb, 0)
        
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef int c_build_tree(self, obbnode_t obb, int level):
        cdef:
            obbnode_t child1, child2
            obbnode_t *parent
            int i, j, split_plane, icell 
            unsigned long parent_idx
            vector_t p, n, c, pt
            double best_ratio, ratio
            int[:,:] cells = self.cells
            int[:] cell
            int[:] cell_ids = self.cell_ids
            int split_acceptable=0, this_cell_id, next_cell_id, start_idx=obb.cell_list_idx, n_cells = obb.n_cells
            int negative, positive, out_cells_left=0, out_cells_right=(n_cells-1)
            
        #Copy the node into the node-list
        parent_idx = self.add_node_c(obb)
        #Get a refrence to the node in the list 
        parent = self.nodes + parent_idx#self.nodes[parent_idx]
        
        if (level < self.max_level) and (n_cells > self.number_of_cells_per_node):
            
            p.x = obb._corner.x + obb._axes[0].x/2. + obb._axes[1].x/2. + obb._axes[2].x/2.
            p.y = obb._corner.y + obb._axes[0].y/2. + obb._axes[1].y/2. + obb._axes[2].y/2.
            p.z = obb._corner.z + obb._axes[0].z/2. + obb._axes[1].z/2. + obb._axes[2].z/2.
                 
            best_ratio = 1.0
            
            for split_plane in range(3):
                n = norm_(obb._axes[split_plane])
                out_cells_left=0
                out_cells_right=0 #out_cell_right + (n_cells-1)
                
                next_cell_id = this_cell_id = cell_ids[start_idx]
                while (out_cells_left + out_cells_right) < n_cells:
                    cell = cells[this_cell_id]
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
                            cell_ids[start_idx + out_cells_left] = this_cell_id
                            out_cells_left += 1
                            next_cell_id = cell_ids[start_idx + out_cells_left] 
                        else:
                            next_cell_id = cell_ids[start_idx + n_cells - 1 - out_cells_right]
                            cell_ids[start_idx + n_cells - 1 - out_cells_right] = this_cell_id
                            out_cells_right += 1
                    else:
                        if negative:
                            cell_ids[start_idx + out_cells_left] = this_cell_id
                            out_cells_left += 1
                            next_cell_id = cell_ids[start_idx + out_cells_left] 
                        else:
                            next_cell_id = cell_ids[start_idx + n_cells - 1 - out_cells_right]
                            cell_ids[start_idx + n_cells - 1 - out_cells_right] = this_cell_id
                            out_cells_right += 1
                    this_cell_id = next_cell_id
                            
                ratio = fabs( <double>(out_cells_right - out_cells_left)/n_cells )
                
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
                child1 = self.compute_obb_cells_c(start_idx, out_cells_left)
                child2 = self.compute_obb_cells_c(start_idx + out_cells_left, out_cells_right)
                child1.parent = parent_idx
                child2.parent = parent_idx
                parent.child1 = self.c_build_tree(child1, level+1)
                parent.child2 = self.c_build_tree(child2, level+1)
        
        return parent_idx 
            
        
        
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
            double[:,:] _axes = self._axes
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
        if _jacobi_no_eigenvals(a, _axes, 1e-9) < 0:
            raise Exception("Jacobi iteration failed.")
        ### Jacobi method outputs unit-magnitude eigenvectors so we don't need to re-normalise
        
        for i in range(3):
            node._axes[i].x = _axes[i,0]
            node._axes[i].y = _axes[i,1]
            node._axes[i].z = _axes[i,2]
        
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
        corner = addvv_(centre, multvs_(node._axes[0], tmin[0]))
        corner = addvv_(corner, multvs_(node._axes[1], tmin[1]))
        corner = addvv_(corner, multvs_(node._axes[2], tmin[2]))
            
        for j in range(3):
            node._axes[j].x *= (tmax[j] - tmin[j])
            node._axes[j].y *= (tmax[j] - tmin[j])
            node._axes[j].z *= (tmax[j] - tmin[j])
                #pt0 = node.axes[i]
                #pt0[j] *= (tmax[j] - tmin[j])
        
        node._corner = corner
        return node
    
    cdef void clear_point_mask(self):
        cdef:
            int i, n=self.point_mask.shape[0]
        for i in range(n):
            self.point_mask[i] = 0
    
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef obbnode_t compute_obb_cells_c(self, int cell_list_idx, int n_cells):
        """
        cell_list - a 1d array of indices into the cells member.
        
        returns:
            - a OBBNode instance
        """
        cdef:
            vector_t mean, p, q, r, dp0, dp1, c
            int i,j, n_pts=self.points.shape[0]
            double tot_mass=0.0, tri_mass
            int[:] tri
            int[:] cell_ids = self.cell_ids
            double[3] _mean
            double[3] tmin = [DBL_MAX, DBL_MAX, DBL_MAX]
            double[3] tmax = [-DBL_MAX, -DBL_MAX, -DBL_MAX]
            ###co-variance matrix
            double[:,:] a = self._covar
            double[:,:] _axes = self._axes
            double[:] a0 = a[0]
            double[:] a1 = a[1]
            double[:] a2 = a[2]
            obbnode_t node
            
        mean.x=0.0
        mean.y=0.0
        mean.z=0.0
        
        node.cell_list_idx = cell_list_idx
        node.n_cells = n_cells
        
        self.clear_point_mask()
            
        for i in range(n_cells):
            tri = self.cells[cell_ids[cell_list_idx + i]]
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
                a[i,j] = (a[i,j] / tot_mass) - (_mean[i] * _mean[j])
                
        ### Find eigenvectors
        if _jacobi_no_eigenvals(a, _axes, 1e-9) < 0:
            raise Exception("Jacobi iteration failed.")
        ### Jacobi method outputs unit-magnitude eigenvectors so we don't need to re-normalise
        
        for i in range(3):
            node._axes[i].x = _axes[i,0]
            node._axes[i].y = _axes[i,1]
            node._axes[i].z = _axes[i,2]
        
        ### Get bounding box by projecting onto eigenvectors
        for i in range(n_pts):
            if self.point_mask[i]:
                vpt = as_vect(self.points[i])
                vpt = subvv_(vpt, mean)
                for j in range(3): #iterate over axes
                    ax = node._axes[j]
                    t = dotprod_(ax, vpt)
                    if t < tmin[j]:
                        tmin[j] = t
                    if t > tmax[j]:
                        tmax[j] = t
        
        #print("min: ", tmin[0], tmin[1], tmin[2], "    max: ", tmax[0], tmax[1], tmax[2])
        corner = addvv_(mean, multvs_(node._axes[0], tmin[0]))
        corner = addvv_(corner, multvs_(node._axes[1], tmin[1]))
        corner = addvv_(corner, multvs_(node._axes[2], tmin[2]))
            
        for j in range(3):
            node._axes[j].x *= (tmax[j] - tmin[j])
            node._axes[j].y *= (tmax[j] - tmin[j])
            node._axes[j].z *= (tmax[j] - tmin[j])
                #pt0 = node.axes[i]
                #pt0[j] *= (tmax[j] - tmin[j])
        
        node._corner = corner
        node.parent = -1
        node.child1 = -1
        node.child2 = -1
        return node
        
            
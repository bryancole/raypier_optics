#cython: language_level=3
"""
Implementation of a Oriented Boundary Boxx Tree (OBB Tree) spatial search algorithm.
Somewhat copied from the vtkOBBTree implementation.
"""
from PIL.TiffImagePlugin import II

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
            addvv_, multvs_, cross_, mag_, Face, intersect_t


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


        
cdef class OBBTree(object):
    def __dealloc__(self):
        free(self.nodes)
        
    def __cinit__(self, double[:,:] points, int[:,:] cells):
        cdef:
            unsigned long max_n_nodes
            
        self.points = points
        self.cells = cells
        self.cell_ids = np.arange(cells.shape[0], dtype=np.int32)
        
        self.point_mask = np.zeros(points.shape[0], dtype=np.int32)
        self.level = 0
        self.max_level = 12
        self.number_of_cells_per_node = 1
        self._covar = np.zeros(shape=(3,3))
        self._axes = np.zeros(shape=(3,3))
        
        max_n_nodes = 2**(int(np.ceil(np.log2(cells.shape[0]))) + 1)
        self.max_n_nodes = max_n_nodes
        self.nodes = <obbnode_t*>malloc(max_n_nodes*sizeof(obbnode_t))
        self.n_nodes = 0
        self.tolerance = 0.1
        
    cdef obbnode_t get_node_c(self, unsigned long i):
        return self.nodes[i]
    
    cdef void set_node_c(self, unsigned long i, obbnode_t node):
        self.nodes[i] = node
        
    cdef unsigned long add_node_c(self, obbnode_t r):       
        cdef:
            unsigned long n_nodes = self.n_nodes     
        if self.n_nodes == self.max_n_nodes:
            if self.max_n_nodes == 0:
                self.max_n_nodes = 1
            else:
                self.max_n_nodes *= 2
            self.nodes = <obbnode_t*>realloc(self.nodes, self.max_n_nodes*sizeof(obbnode_t))
        self.nodes[n_nodes] = r
        self.n_nodes += 1
        return n_nodes
    
    property root:
        def __get__(self):
            cdef:
                long root_idx = self.root_idx
            if root_idx<0:
                return None
            root = OBBNode()
            root.owner = self
            root.node = self.get_node_c(root_idx)
            return root
        
    def clear_tree(self):
        self.n_nodes = 0
        self.level = 0
        
    def build_tree(self):
        cdef:
            int n_cells = self.cells.shape[0]
            obbnode_t obb
            long root_idx
            OBBNode root

        self.n_nodes = 0
        self.level = 0
        obb = self.compute_obb_cells_c(0, n_cells)
        self.root_idx = self.c_build_tree(obb, 0)
        
    cdef int line_intersects_node_c(self, obbnode_t* this, vector_t p0, vector_t p1):
        cdef:
            int ii
            #obbnode_t this = self.nodes[node_idx]
            double rangeAmin, rangeAmax, rangeBmin, rangeBmax
            double eps = self.tolerance
            
        for ii in range(3):
            rangeAmin = dotprod_(this._corner, this._axes[ii])
            rangeAmax = rangeAmin + dotprod_(this._axes[ii], this._axes[ii])
            
            rangeBmin = dotprod_(p0, this._axes[ii])
            rangeBmax = rangeBmin
            dotB = dotprod_(p1, this._axes[ii])
            if dotB < rangeBmin:
                rangeBmin = dotB
            else:
                rangeBmax = dotB
                
            if eps != 0:
                eps *= sqrt(fabs(rangeAmax - rangeAmin))
                
            if ((rangeAmax + eps < rangeBmin) | (rangeBmax + eps < rangeAmin)):
                return 0
        return 1
    
    @boundscheck(False)
    @initializedcheck(False)
    cdef vector_t get_point_c(self, long idx):
        cdef:
            vector_t p
        p.x = self.points[idx,0]
        p.y = self.points[idx,1]
        p.z = self.points[idx,2]
        return p
    
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef double line_intersects_cell_c(self, long cell_idx, vector_t o, vector_t d):
        """ For a line defined by origin 'o' and direction 'd',
        returns the distance form the origin to the triangle intersection.
        If there is no intersection, returns -1
        """
        cdef:
            int[:] cell=self.cells[cell_idx,:]
            vector_t n, p1, p2, p3
            vector_t v1, v2, a0, da0
            double alpha, det, invdet, u,v
            
        p1 = self.get_point_c(cell[0])
        p2 = self.get_point_c(cell[1])
        p3 = self.get_point_c(cell[2])
            
        v1 = subvv_(p2,p1)
        v2 = subvv_(p3,p1)
        n = cross_(v1, v2)
        
        det = -dotprod_(d,n)
        if det == 0.0:
            return -1
            #print(d.x,d.y,d.z,n.x,n.y,n.z)
        invdet = 1.0/det
        
        a0 = subvv_(o, p1)
        da0 = cross_(a0, d)
        u = dotprod_(v2, da0) * invdet
        v = -dotprod_(v1, da0) * invdet
        alpha = dotprod_(a0,n) * invdet
        
        if (u+v > 1.0) | (u<0) | (v<0) | (alpha<0) :
            return -1.0
        
        return alpha
        
        
#         bool intersect_triangle(
#     in Ray R, in vec3 A, in vec3 B, in vec3 C, out float t, 
#     out float u, out float v, out vec3 N
# ) { 
#    vec3 E1 = B-A;
#    vec3 E2 = C-A;
#          N = cross(E1,E2);
#    float det = -dot(R.Dir, N);
#    float invdet = 1.0/det;
#    vec3 AO  = R.Origin - A;
#    vec3 DAO = cross(AO, R.Dir);
#    u =  dot(E2,DAO) * invdet;
#    v = -dot(E1,DAO) * invdet;
#    t =  dot(AO,N)  * invdet; 
#    return (det >= 1e-6 && t >= 0.0 && u >= 0.0 && v >= 0.0 && (u+v) <= 1.0);
# }

    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef intersection intersect_with_line_c(self, vector_t p1, vector_t p2, long long[:] workspace):
        """
        """
        cdef:
            long depth=1, cell_idx
            obbnode_t *node
            vector_t d = subvv_(p2,p1)
            double alpha
            intersection res
            int i
            double tol=self.tolerance/mag_(d)
            
        res.alpha = 1.0
        res.cell_idx = -1
        workspace[0] = self.root_idx
        while depth > 0:
            depth -= 1
            node = self.nodes + workspace[depth]
            if self.line_intersects_node_c(node, p1, p2):
                if (node.child1 < 0) or (node.child2 < 0):
                    for i in range(node.n_cells):
                        cell_idx = self.cell_ids[node.cell_list_idx+i]
                        alpha = self.line_intersects_cell_c(cell_idx, p1, d)
                        if (alpha >= tol) and (alpha < res.alpha):
                            res.alpha = alpha
                            res.cell_idx = cell_idx
                else:
                    workspace[depth] = node.child1
                    workspace[depth+1] = node.child2
                    depth += 2
        #print("cell intersection", icount, "obb intersections", obb_count, "total", icount+obb_count)
        if res.cell_idx < 0:
            res.alpha = -1.0
        return res 
        
        
    def intersect_with_line(self, point1, point2):
        """
        Find the closest intersection to point 1 with the line specified by the
        given points.
        """
        cdef:
            vector_t o,d
            intersection its
            double alpha
            vector_t ip
            
        o.x = point1[0]
        o.y = point1[1]
        o.z = point1[2]
        d.x = point2[0]
        d.y = point2[1]
        d.z = point2[2]
        if self.level < 1:
            raise ValueError("OBBTree structure not yet built.")
        workspace = np.zeros(self.level+1, dtype=np.int64)
        its = self.intersect_with_line_c(o,d, workspace)
        if its.cell_idx < 0:
            return (None, -1)
        else:
            ip = addvv_(o, multvs_(subvv_(d,o), its.alpha))
            return ((ip.x, ip.y, ip.z), its.cell_idx)
        
        
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef long c_build_tree(self, obbnode_t obb, int level):
        cdef:
            obbnode_t child1, child2
            obbnode_t *parent
            int i, j, split_plane, icell 
            long parent_idx
            vector_t p, n, c, pt
            double best_ratio, ratio
            int[:,:] cells = self.cells
            int[:] cell
            int[:] cell_ids = self.cell_ids
            int split_acceptable=0, found_best_split=0 
            int this_cell_id, next_cell_id, start_idx=obb.cell_list_idx, n_cells = obb.n_cells
            int negative, positive, out_cells_left=0, out_cells_right=(n_cells-1)
            
        #Copy the node into the node-list
        parent_idx = self.add_node_c(obb)
        #Get a refrence to the node in the list 
        parent = self.nodes + parent_idx#self.nodes[parent_idx]
        
        if level > self.level:
            self.level = level
        
        if (level < self.max_level) and (n_cells > self.number_of_cells_per_node):
            
            p.x = obb._corner.x + obb._axes[0].x/2. + obb._axes[1].x/2. + obb._axes[2].x/2.
            p.y = obb._corner.y + obb._axes[0].y/2. + obb._axes[1].y/2. + obb._axes[2].y/2.
            p.z = obb._corner.z + obb._axes[0].z/2. + obb._axes[1].z/2. + obb._axes[2].z/2.
                 
            best_ratio = 1.0
            
            split_plane = 0
            while (split_plane < 3) and not split_acceptable:
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
                            
                ratio = fabs( (<double>(out_cells_right - out_cells_left))/n_cells )
                
                if (ratio <= 0.6) or found_best_split:
                    split_acceptable = 1
                    break
                else:
                    if ratio < best_ratio:
                        best_ratio = ratio
                        best_plane = split_plane
                        
                    split_plane += 1
                    
                    if (split_plane==3) and (best_ratio < 0.95):
                        split_plane = best_plane
                        found_best_split = 1
                        
                
            if split_acceptable:
                #print(level, ratio, n_cells, out_cells_left, out_cells_right)
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
        
        ## reset covariance matrix
        for i in range(3):
            a0[i] = 0.0
            a1[i] = 0.0
            a2[i] = 0.0
            
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
        
        
cdef class OBBNode(object):   
    property parent:
        def __get__(self):
            cdef:
                long parent = self.node.parent
            if parent<0:
                return None
            node = OBBNode()
            node.owner = self.owner
            node.node = self.owner.get_node_c(parent)
            return node
        
    property child1:
        def __get__(self):
            cdef:
                long child = self.node.child1
            if child<0:
                return None
            node = OBBNode()
            node.owner = self.owner
            node.node = self.owner.get_node_c(child)
            return node
        
    property child2:
        def __get__(self):
            cdef:
                long child = self.node.child2
            if child<0:
                return None
            node = OBBNode()
            node.owner = self.owner
            node.node = self.owner.get_node_c(child)
            return node
        
    property corner:
        def __get__(self):
            return np.array([self.node._corner.x, self.node._corner.y, self.node._corner.z])
        
        def __set__(self, c):
            self.node._corner.x = c[0]
            self.node._corner.y = c[1]
            self.node._corner.z = c[2]
            
    property axes:
        def __get__(self):
            cdef:
                int i
            out = np.empty(shape=(3,3))
            for i in range(3):
                out[i,0] = self.node._axes[i].x
                out[i,1] = self.node._axes[i].y
                out[i,2] = self.node._axes[i].z
            return out
        
        def __set__(self, double[:,:] axes):
            cdef:
                int i
            for i in range(3):
                self.node._axes[i].x = axes[i,0]
                self.node._axes[i].y = axes[i,1]
                self.node._axes[i].z = axes[i,2]
                    
    @boundscheck(False)
    @initializedcheck(False)
    @cdivision(True)
    cdef int intersect_line_c(self, vector_t p1, vector_t p2):
        cdef:
            int i
            double rangeAmin, rangeAmax, rangePmin, rangePmax, dotP
            vector_t ax
            
        for i in range(3):
            ax = self.node._axes[i]
            
            rangeAmin = dotprod_(self.node._corner, ax)
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
                    
                
cdef class OBBTreeFace(Face):
    cdef:
        public OBBTree obbtree
        long long[:] workspace
        double[:,:] normals
        
    def __cinit__(self, **kwds):
        cdef:
            OBBTree tree = kwds['tree']
            double[:,:] normals
            int i, n_cells = tree.cells.shape[0]
            vector_t n, p1, p2, p3
            int[:] cell
            
        self.obbtree = tree
        if tree.level <= 0:
            tree.build_tree()
        self.workspace = np.zeros(tree.level+1, np.int64)
        normals = np.empty((n_cells,3), np.float64)
        
        for i in range(n_cells):
            cell = tree.cells[i]
            p1 = tree.get_point_c(cell[0])
            p2 = tree.get_point_c(cell[1])
            p3 = tree.get_point_c(cell[2])
            n = norm_(cross_(subvv_(p2,p1), subvv_(p3,p1)))
            normals[i,0] = n.x
            normals[i,1] = n.y
            normals[i,2] = n.z
            
        self.normals = normals 
        
        
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
            intersection it
            intersect_t out
            
        it = self.obbtree.intersect_with_line_c(p1, p2, self.workspace)
        
        out.dist = it.alpha * mag_(subvv_(p2,p1))
        out.piece_idx = <int>(it.cell_idx) #casting long to int. Hopefully there are not too many cells
        return out
        
        
    cdef vector_t compute_normal_c(self, vector_t p, int piece):
        """Compute the surface normal in local coordinates,
        given a point on the surface (also in local coords).
        """
        cdef:
            vector_t n
            double[:] na
        
        na = self.normals[piece]
        n.x = na[0]
        n.y = na[1]
        n.z = na[2]
        return n
    
    
                
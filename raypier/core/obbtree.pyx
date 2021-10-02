
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
    
import numpy as np
cimport numpy as np_


cdef struct MaxElem:
    double Amax
    int i
    int j 

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


def jacobi(double[:,:] A, double tol=1e-9):
    cdef:
        int i, n = len(A)
        int max_rot = 5*(n**2)
        double[:,:] p = np.identity(n)*1.0
        MaxElem me
        
    for i in range(max_rot):
        me = max_elem(A)
        if me.Amax < tol:
            return np.diagonal(A),p
        rotate(A,p,me.i, me.j)
    raise Exception("Jacobi method did not converge")



cdef class OBBNode(object):
    cdef:
        double[3] centre #Centre coordinate of the box
        double[3][3] axes #the three primary axes of the box, ordered from longest to smallest
        OBBNode parent, child1, child2 #children will be None for leaf nodes. Parent will be None for root
        int[:] cell_list
        
        
        
cdef class OBBTree(object):
    cdef:
        double[:,:] points
        int[:,:] cells #always triangles
        OBBNode root
        
    def __init__(self, double[:,:] points, int[:,:] cells):
        self.points = points
        self.cells = cells
        
    def compute_obb(self, int[:] point_ids):
        """
        Compute an OBB from the list of points given. Return the corner point
        and the three axes defining the orientation of the OBB. Also return
        a sorted list of relative "sizes" of axes for comparison purposes.
        
        In fact, we want to compute the OBB based on cell 'moments'. These are 
        weighed by cell area.
        """
        cdef:
            double[3] centre=[0,0,0], pt0
            double[:] pt
            int i,j, n_pts = point_ids.size
            ###co-variance matrix
            double[3][3] a = [[0.0,0.0,0.0],
                              [0.0,0.0,0.0],
                              [0.0,0.0,0.0]]
        
        for i in range(n_pts):
            pt = self.points[point_ids[i]]
            for j in range(3):
                centre[j] += pt[j]
        for j in range(3):
            centre[j] /= n_pts

        ### Populate covariance matrix
        for i in range(n_pts):
            pt = self.points[point_ids[i]]
            pt0[0] = pt[0] - centre[0]
            pt0[1] = pt[1] - centre[1]
            pt0[2] = pt[2] - centre[2]
            for j in range(3):
                a[0][j] += pt0[0] * pt0[j]
                a[1][j] += pt0[1] * pt0[j]
                a[2][j] += pt0[2] * pt0[j]
        for j in range(3):
            a[0][j] /= n_pts
            a[1][j] /= n_pts
            a[2][j] /= n_pts
            
        ### Find eigenvectors
        
        


from .ctracer cimport vector_t, subvv_, dotprod_, mag_sq_, norm_,\
            addvv_, multvs_, cross_, mag_, Face, intersect_t

cdef struct intersection:
    double alpha
    long cell_idx


cdef struct obbnode_t:
    vector_t _corner #one corner of the box
    vector_t[3] _axes #the three primary axes of the box, ordered from longest to smallest
    long parent, child1, child2 #children will be -1 for leaf nodes. Parent will be -1 for root
    int cell_list_idx
    int n_cells


cdef class OBBTree(object):
    cdef:
        double[:,:] points
        int[:,:] cells #always triangles
        int[:] cell_ids
        int[:] point_mask
        
        long root_idx
        
        public int level
        public int max_level
        public int number_of_cells_per_node
        
        double[:,:] _covar #A workspace for obb calculation
        double[:,:] _axes
        
        unsigned long n_nodes
        unsigned long max_n_nodes
        obbnode_t *nodes
        
        public double tolerance
        
    cdef obbnode_t get_node_c(self, unsigned long i)
    cdef void set_node_c(self, unsigned long i, obbnode_t node)
    cdef unsigned long add_node_c(self, obbnode_t r)
    cdef int line_intersects_node_c(self, obbnode_t* this, vector_t p0, vector_t p1)
    cdef vector_t get_point_c(self, long idx)
    cdef double line_intersects_cell_c(self, long cell_idx, vector_t o, vector_t d)
    cdef intersection intersect_with_line_c(self, vector_t p1, vector_t p2, long long[:] workspace)
    cdef long c_build_tree(self, obbnode_t obb, int level)
    cdef void clear_point_mask(self)
    cdef obbnode_t compute_obb_cells_c(self, int cell_list_idx, int n_cells)


cdef class OBBNode(object):
    cdef:
        OBBTree owner
        obbnode_t node
        
    cdef int intersect_line_c(self, vector_t p1, vector_t p2)
    
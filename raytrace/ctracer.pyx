#!/bin/env python

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

cdef extern from "math.h":
    double sqrt(double arg)
    double fabs(double arg)
    #double INFINITY
    
cdef extern from "float.h":
    double DBL_MAX

cdef INF=(DBL_MAX+DBL_MAX)
    
from libc.stdlib cimport malloc, free, realloc

cdef extern from "stdlib.h" nogil:
    void *memcpy(void *str1, void *str2, size_t n)

import numpy as np
cimport numpy as np_

ray_dtype = np.dtype([('origin', np.double, (3,)),
                        ('direction', np.double, (3,)),
                        ('normal', np.double, (3,)),
                        ('E_vector', np.double, (3,)),
                        ('refractive_index', np.complex128),
                        ('E1_amp', np.complex128),
                        ('E2_amp', np.complex128),
                        ('length', np.double),
                        ('wavelength', np.double),
                        ('parent_idx', np.uint32),
                        ('end_face_idx', np.uint32)
                        ])
                        

############################################
### C type declarations for internal use ###
############################################

cdef struct vector_t:
    double x,y,z

    
cdef struct complex_t:
    double real
    double imag

    
cdef struct ray_t:
    #vectors
    vector_t origin, direction, normals, E_vector
    #complex attribs
    complex_t refractive_index, E1_amp, E2_amp
    #simple attribs
    double length, wavelength
    #reference ids to related objects
    unsigned int parent_idx, end_face_idx
    ##objects
    #object face, end_face, child_refl, child_trans
    
cdef struct transform_t:
    double m00, m01, m02, m10, m11, m12, m20, m21, m22
    double tx, ty, tz


##############################
### Vector maths functions ###
##############################

cdef inline vector_t transform_c(transform_t t, vector_t p):
    cdef vector_t out
    out.x = p.x*t.m00 + p.y*t.m01 + p.z*t.m02 + t.tx
    out.y = p.x*t.m10 + p.y*t.m11 + p.z*t.m12 + t.ty
    out.z = p.x*t.m20 + p.y*t.m21 + p.z*t.m22 + t.tz
    return out

cdef inline vector_t rotate_c(transform_t t, vector_t p):
    cdef vector_t out
    out.x = p.x*t.m00 + p.y*t.m01 + p.z*t.m02
    out.y = p.x*t.m10 + p.y*t.m11 + p.z*t.m12
    out.z = p.x*t.m20 + p.y*t.m21 + p.z*t.m22
    return out

cdef inline vector_t set_v(object O):
    cdef vector_t v
    v.x = O[0]
    v.y = O[1]
    v.z = O[2]
    return v

def py_set_v(O):
    cdef vector_t v_
    v_ = set_v(O)
    return (v_.x, v_.y, v_.z)

cdef inline double sep_(vector_t p1, vector_t p2):
    cdef double a,b
    a = (p2.x-p1.x)
    b = (p2.y-p1.y)
    c = (p2.z-p1.z)
    return sqrt((a*a) + (b*b) + (c*c))

def sep(a, b):
    cdef vector_t a_ = set_v(a), b_ = set_v(b)
    return sep_(a_, b_)

cdef inline vector_t invert_(vector_t v):
    v.x = -v.x
    v.y = -v.y
    v.z = -v.z
    return v

def invert(v):
    cdef vector_t v_ = set_v(v)
    v_ = invert_(v_)
    return (v_.x, v_.y, v_.z)

cdef inline vector_t multvv_(vector_t a, vector_t b):
    cdef vector_t out
    out.x = a.x*b.x
    out.y = a.y*b.y
    out.z = a.z*b.z
    return out

def multvv(a, b):
    cdef vector_t a_, b_, c_
    a_ = set_v(a)
    b_ = set_v(b)
    c_ = multvv_(a_, b_)
    return (c_.x, c_.y, c_.z)

cdef inline vector_t multvs_(vector_t a, double b):
    cdef vector_t out
    out.x = a.x*b
    out.y = a.y*b
    out.z = a.z*b
    return out

def multvs(a, b):
    cdef vector_t a_, c_
    a_ = set_v(a)
    c_ = multvs_(a_, b)
    return (c_.x, c_.y, c_.z)

cdef inline vector_t addvv_(vector_t a, vector_t b):
    cdef vector_t out
    out.x = a.x+b.x
    out.y = a.y+b.y
    out.z = a.z+b.z
    return out

def addvv(a, b):
    cdef vector_t a_, b_, c_
    a_ = set_v(a)
    b_ = set_v(b)
    c_ = addvv_(a_, b_)
    return (c_.x, c_.y, c_.z)
    
cdef inline vector_t addvs_(vector_t a, double b):
    cdef vector_t out
    out.x = a.x+b
    out.y = a.y+b
    out.z = a.z+b
    return out

def addvs(a, b):
    cdef vector_t a_, c_
    a_ = set_v(a)
    c_ = addvs_(a_, b)
    return (c_.x, c_.y, c_.z)

cdef inline vector_t subvv_(vector_t a, vector_t b):
    cdef vector_t out
    out.x = a.x-b.x
    out.y = a.y-b.y
    out.z = a.z-b.z
    return out

def subvv(a, b):
    cdef vector_t a_, b_, c_
    a_ = set_v(a)
    b_ = set_v(b)
    c_ = subvv_(a_, b_)
    return (c_.x, c_.y, c_.z)

cdef inline vector_t subvs_(vector_t a, double b):
    cdef vector_t out
    out.x = a.x-b
    out.y = a.y-b
    out.z = a.z-b
    return out

def subvs(a, b):
    cdef vector_t a_, c_
    a_ = set_v(a)
    c_ = subvs_(a_, b)
    return (c_.x, c_.y, c_.z)

cdef inline double mag_(vector_t a):
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z)

def mag(a):
    cdef vector_t a_
    a_ = set_v(a)
    return mag_(a_)

cdef inline double mag_sq_(vector_t a):
    return a.x*a.x + a.y*a.y + a.z*a.z

def mag_sq(a):
    cdef vector_t a_
    a_ = set_v(a)
    return mag_sq_(a_)

cdef inline double dotprod_(vector_t a, vector_t b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def dotprod(a, b):
    cdef vector_t a_, b_
    a_ = set_v(a)
    b_ = set_v(b)
    return dotprod_(a_,b_)

cdef inline vector_t cross_(vector_t a, vector_t b):
    cdef vector_t c
    c.x = a.y*b.z - a.z*b.y
    c.y = a.z*b.x - a.x*b.z
    c.z = a.x*b.y - a.y*b.x
    return c

def cross(a, b):
    cdef vector_t a_, b_, c_
    a_ = set_v(a)
    b_ = set_v(b)
    c_ = cross_(a_, b_)
    return (c_.x, c_.y, c_.z)

cdef vector_t norm_(vector_t a):
    cdef double mag=sqrt(a.x*a.x + a.y*a.y + a.z*a.z)
    a.x /= mag
    a.y /= mag
    a.z /= mag
    return a

def norm(a):
    cdef vector_t a_
    a_ = set_v(a)
    a_ = norm_(a_)
    return (a_.x, a_.y, a_.z)

##################################
### Python extension types
##################################

cdef class Transform:
    
    def __init__(self, rotation=[[1,0,0],[0,1,0],[0,0,1]], 
                        translation=[0,0,0]):
        self.rotation = rotation
        self.translation = translation
        
    property rotation:
        def __set__(self, rot):
            cdef transform_t t
            t.m00, t.m01, t.m02 = rot[0]
            t.m10, t.m11, t.m12 = rot[1]
            t.m20, t.m21, t.m22 = rot[2]
            self.trans = t
            
        def __get__(self):
            cdef transform_t t
            t = self.trans
            return [[t.m00, t.m01, t.m02],
                    [t.m10, t.m11, t.m12],
                    [t.m20, t.m21, t.m22]]
    
    property translation:
        def __set__(self, dt):
            self.trans.tx, self.trans.ty, self.trans.tz = dt
        
        def __get__(self):
            return (self.trans.tx, self.trans.ty, self.trans.tz)


cdef class RayCollectionIterator:        
    def __cinit__(self, RayCollection rays):
        self.rays = rays
        self.counter = 0
        
    def __iter__(self):
        return self
    
    def __next__(self):
        cdef Ray ray=Ray.__new__(Ray)
        if self.counter >= self.rays.n_rays:
            raise StopIteration
        ray.ray = self.rays.rays[self.counter]
        self.counter += 1
        return ray


cdef class Ray:
    
    def __cinit__(self, **kwds):
        for k in kwds:
            setattr(self, k, kwds[k])
            
    def __repr__(self):
        return "Ray(o=%s, d=%s)"%(str(self.origin),
                                            str(self.direction))
                
    property origin:
        """Origin coordinates of the ray"""
        def __get__(self):
            return (self.ray.origin.x,self.ray.origin.y,self.ray.origin.z)
        
        def __set__(self, v):
            self.ray.origin.x = v[0]
            self.ray.origin.y = v[1]
            self.ray.origin.z = v[2]
            
    property direction:
        """direction of the ray, normalised to a unit vector"""
        def __get__(self):
            return (self.ray.direction.x,self.ray.direction.y,self.ray.direction.z)
        
        def __set__(self, v):
            self.ray.direction.x = v[0]
            self.ray.direction.y = v[1]
            self.ray.direction.z = v[2]
            
    property normals:
        """normal vector for the face which created this ray"""
        def __get__(self):
            return (self.ray.normal.x,self.ray.normal.y,self.ray.normal.z)
        
        def __set__(self, v):
            self.ray.normal.x = v[0]
            self.ray.normal.y = v[1]
            self.ray.normal.z = v[2]
            
    property E_vector:
        """Unit vector, perpendicular to the ray direction,
        which gives the direction of E-field polarisation"""
        def __get__(self):
            return (self.ray.E_vector.x,self.ray.E_vector.y,self.ray.E_vector.z)
        
        def __set__(self, v):
            self.ray.E_vector.x = v[0]
            self.ray.E_vector.y = v[1]
            self.ray.E_vector.z = v[2]
            
    property length:
        """The length of the ray. This is infinite in 
        unterminated rays"""
        def __get__(self):
            return self.ray.length
        
        def __set__(self, double v):
            self.ray.length = v
            
    property termination:
        """the end-point of the ray (read only)
        """
        def __get__(self):
            cdef vector_t end
            cdef float length
            if self.ray.length > 100.0:
                length = 100.0
            else:
                length = self.ray.length
            end = addvv_(self.ray.origin, multvs_(self.ray.direction, 
                                    length))
            return (end.x, end.y, end.z)
        
    property refractive_index:
        """complex refractive index through which
        this ray is propagating"""
        def __get__(self):
            return complex(self.ray.refractive_index.real,
                            self.ray.refractive_index.imag)
        def __set__(self, v):
            self.ray.refractive_index.real = v.real
            self.ray.refractive_index.imag = v.imag
            
    property E1_amp:
        """Complex amplitude of the electric field polarised
        parallel to the E_vection."""
        def __get__(self):
            return complex(self.ray.E1_amp.real,
                            self.ray.E1_amp.imag)
        def __set__(self, v):
            self.ray.E1_amp.real = v.real
            self.ray.E1_amp.imag = v.imag
            
    property E2_amp:
        """Complex amplitude of the electric field polarised
        perpendicular to the E_vection"""
        def __get__(self):
            return complex(self.ray.E2_amp.real,
                            self.ray.E2_amp.imag)
        def __set__(self, v):
            self.ray.E2_amp.real = v.real
            self.ray.E2_amp.imag = v.imag
            
    property parent_idx:
        """Index of the parent ray in parent RayCollection
        """
        def __get__(self):
            return self.ray.parent_idx
        
        def __set__(self, int v):
            self.ray.parent_idx = v
            
    property end_face_idx:
        """Index of the terminating face, in the global
        face list (created for each tracing operation)
        """
        def __get__(self):
            return self.ray.end_face_idx
        
        def __set__(self, int v):
            self.ray.end_face_idx = v


cdef class RayCollection:
    
    def __cinit__(self, size_t max_size):
        self.rays = <ray_t*>malloc(max_size*sizeof(ray_t))
        self.n_rays = 0
        self.max_size = max_size
        
    def __dealloc__(self):
        free(self.rays)
        
    def __len__(self):
        return self.n_rays
        
    cdef add_ray_c(self, ray_t r):
        if self.n_rays == self.max_size:
            if self.max_size == 0:
                self.max_size = 1
            else:
                self.max_size *= 2
            self.rays = <ray_t*>realloc(self.rays, self.max_size*sizeof(ray_t))
        self.rays[self.n_rays] = r
        self.n_rays += 1
        
    def reset_length(self):
        """Sets the length of all rays in this RayCollection to Infinity
        """
        cdef int i
        for i in xrange(self.n_rays):
            self.rays[i].length = INF
        
    def add_ray(self, Ray r):
        """Adds the given Ray instance to this collection
        """
        self.add_ray_c(r.ray)
        
    def add_ray_list(self, list rays):
        """Adds the given list of Rays to this collection
        """
        cdef int i
        for i in xrange(len(rays)):
            if not isinstance(rays[i], Ray):
                raise TypeError("ray list contains non-Ray instance at index %d"%i)
        for i in xrange(len(rays)):
            self.add_ray_c((<Ray>rays[i]).ray)
        
    def clear_ray_list(self):
        """Empties this RayCollection (by setting the count to zero)
        """
        self.n_rays = 0
        
    def get_ray_list(self):
        """Returns the contents of this RayCollection as a list of Rays
        """
        cdef int i
        cdef list ray_list = []
        cdef Ray r
        for i in xrange(self.n_rays):
            r = Ray()
            r.ray = self.rays[i]
            ray_list.append(r)
        return ray_list
    
    def __getitem__(self, int idx):
        cdef Ray r
        if idx >= self.n_rays:
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        r = Ray()
        r.ray = self.rays[idx]
        return r
    
    def __setitem__(self, int idx, Ray r):
        if idx >= self.n_rays:
            raise IndexError("Attempting to set index %d from a size %d array"%(idx, self.n_rays))
        self.rays[idx] = r.ray
        
    def __iter__(self):
        return RayCollectionIterator(self)
    
    def copy_as_array(self):
        """Returns the contents of this RayCollection as a numpy array
        (the data is always copied).
        """
        cdef np_.ndarray out = np.empty(self.n_rays, dtype=ray_dtype)
        memcpy(<np_.float64_t *>out.data, self.rays, self.n_rays*sizeof(ray_t))
        return out
    
    @classmethod
    def from_array(cls, np_.ndarray data):
        """Creates a new RayCollection from the given numpy array. The array
        dtype should be a ctracer.ray_dtype. The data is copied into the 
        RayCollection
        """
        cdef int size=data.shape[0]
        cdef RayCollection rc = RayCollection(size)
        assert data.dtype is ray_dtype
        memcpy(rc.rays, <np_.float64_t *>data.data, size*sizeof(ray_t))
        rc.n_rays = size
        return rc
        
    
cdef class InterfaceMaterial(object):
    """Abstract base class for objects describing
    the materials characterics of a Face
    """
    
    cdef eval_child_ray_c(self, ray_t *old_ray, 
                                unsigned int ray_idx, 
                                vector_t p, vector_t normal,
                                RayCollection new_rays):
        pass
    
    def eval_child_ray(self, Ray old_ray, ray_idx, point, 
                        normal, RayCollection new_rays):
        cdef:
            vector_t p, n
            Ray out=Ray()
            unsigned int idx
        
        p = set_v(point)
        n = set_v(normal)
        self.eval_child_ray_c(&old_ray.ray, ray_idx, 
                                        p, n, new_rays)
    
    
cdef class Face(object):
    
    params = []
    
    def __cinit__(self, owner=None, tolerance=0.0001, 
                        max_length=100, material=None, **kwds):
        self.name = "base Face class"
        self.tolerance = tolerance
        self.owner = owner
        self.max_length = max_length
        if isinstance(material, InterfaceMaterial):
            self.material = material
        else:
            from raytrace.cmaterials import PECMaterial
            self.material = PECMaterial()
        self.invert_normal = int(kwds.get('invert_normal', 0))
        
    
    cdef double intersect_c(self, vector_t p1, vector_t p2):
        """returns the distance of the nearest valid intersection between 
        p1 and p2. p1 and p2 are in the local coordinate system
        """
        return 0
    
    def update(self):
        """Called to update the parameters from the owner
        to the Face
        """
        for name in self.params:
            v = getattr(self.owner, name)
            setattr(self, name, v)
    
    def intersect(self, p1, p2):
        cdef:
            vector_t p1_, p2_
            double dist
        
        p1_ = set_v(p1)
        p2_ = set_v(p2)
        dist = self.intersect_c(p1_, p2_)
        return dist

    cdef vector_t compute_normal_c(self, vector_t p):
        return p
    
    def compute_normal(self, p):
        """Compute normal vector at a given point, in local
        face coordinates
        """
        cdef vector_t p_, n
        n = self.compute_normal_c(p_)
        return (n.x, n.y, n.z)
        


cdef class FaceList(object):
    """A group of faces which share a transform"""
    def __cinit__(self, owner=None):
        self.transform = Transform()
        self.inverse_transform = Transform()
        self.owner = owner
        
    def sync_transforms(self):
        """sets the transforms from the owner's VTKTransform
        """
        try:
            trans = self.owner.transform
        except AttributeError:
            print "NO OWNER", self.owner
            return
        m = trans.matrix
        rot = [[m.get_element(i,j) for j in xrange(3)] for i in xrange(3)]
        dt = [m.get_element(i,3) for i in xrange(3)]
        #print "TRANS", rot, dt
        self.transform = Transform(rotation=rot, translation=dt)
        inv_trans = trans.linear_inverse
        m = inv_trans.matrix
        rot = [[m.get_element(i,j) for j in xrange(3)] for i in xrange(3)]
        dt = [m.get_element(i,3) for i in xrange(3)]
        self.inverse_transform = Transform(rotation=rot, translation=dt)
        
    property transform:
        def __set__(self, Transform t):
            self.trans = t.trans
           
        def __get__(self):
            cdef Transform t=Transform()
            t.trans = self.trans
            return t
        
    property inverse_transform:
        def __set__(self, Transform t):
            self.inv_trans = t.trans
           
        def __get__(self):
            cdef Transform t=Transform()
            t.trans = self.inv_trans
            return t
        
    def __getitem__(self, intidx):
        return self.faces[intidx]
        
     
    cdef int intersect_c(self, ray_t *ray, vector_t ray_end, double max_length):
        """Finds the face with the nearest intersection
        point, for the ray defined by the two input points,
        P1 and P2 (in global coords).
        """
        cdef:
            vector_t p1 = transform_c(self.inv_trans, ray.origin)
            vector_t p2 = transform_c(self.inv_trans, ray_end)
            list faces=self.faces
            unsigned int i
            int all_idx=-1
            double dist
            Face face
        
        for i in xrange(len(faces)):
            face = faces[i]
            dist = face.intersect_c(p1, p2)
            if face.tolerance < dist < ray.length:
                ray.length = dist
                all_idx = face.idx
                ray.end_face_idx = all_idx
        return all_idx
    
    def intersect(self, Ray r, double max_length):
        cdef vector_t P1_
        cdef int idx
        
        P1_ = addvv_(r.ray.origin, multvs_(r.ray.direction, r.ray.length))
        idx = self.intersect_c(&r.ray, P1_, max_length)
        return idx
    
    cdef vector_t compute_normal_c(self, Face face, vector_t point):
        cdef vector_t out
        
        out = transform_c(self.inv_trans, point)
        out = face.compute_normal_c(out)
        if face.invert_normal:
            out = invert_(out)
        out = rotate_c(self.trans, out)
        return out
    
    def compute_normal(self, Face face, point):
        cdef vector_t p
        p = set_v(point)
        p = self.compute_normal_c(face, p)
        return (p.x, p.y, p.z)
    

##################################
### Python module functions
##################################

cdef RayCollection trace_segment_c(RayCollection rays, 
                                    list face_sets, 
                                    list all_faces,
                                    float max_length):
    cdef:
        FaceList face_set #a FaceList
        unsigned int size, i, j
        vector_t P1, normal, point
        int idx, nearest_set=-1, nearest_idx=-1, n_sets=len(face_sets)
        ray_t new_ray
        ray_t *ray
        RayCollection new_rays
   
    #need to allocate the output rays here 
    new_rays = RayCollection(rays.n_rays)
    
    for i in range(rays.n_rays):
        ray = rays.rays + i
        ray.length = max_length
        ray.end_face_idx = -1
        nearest_idx=-1
        point = addvv_(ray.origin, 
                            multvs_(ray.direction, 
                                    max_length))
        #print "points", P1, P2
        for j in xrange(n_sets):
            face_set = face_sets[j]
            #intersect_c returns the face idx of the intersection, or -1 otherwise
            idx = (<FaceList>face_set).intersect_c(ray, point, max_length)
            if idx >= 0:
                nearest_set = j
                nearest_idx = idx
        if nearest_idx >= 0:
            #print "GET FACE", nearest.face_idx, len(all_faces)
            face = all_faces[nearest_idx]
            face.count += 1
            #print "ray length", ray.length
            point = addvv_(ray.origin, multvs_(ray.direction, ray.length))
            normal = (<FaceList>(face_sets[nearest_set])).compute_normal_c(face, point)
            #print "s normal", normal
            (<InterfaceMaterial>(face.material)).eval_child_ray_c(ray, i, 
                                                    point,
                                                    normal,
                                                    new_rays
                                                    )
    return new_rays


def trace_segment(RayCollection rays, 
                    list face_sets, 
                    list all_faces,
                    max_length=100):
    for fs in face_sets:
        fs.sync_transforms()
    return trace_segment_c(rays, face_sets, all_faces, max_length)


def transform(Transform t, p):
    cdef vector_t p1, p2
    assert isinstance(t, Transform)
    assert len(p)==3
    p1.x = p[0]
    p1.y = p[1]
    p1.z = p[2]
    p2 = transform_c(t.trans, p1)
    return (p2.x, p2.y, p2.z)
    
    

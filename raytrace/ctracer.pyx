#!/usr/bin/env python

cdef extern from "math.h":
    double sqrt(double arg)
    double INFINITY
    
from stdlib cimport malloc, free

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
    
    
cdef struct intersection_t:
    unsigned int face_idx #the intersecting face
    vector_t point #position of intersection
    double dist #fractional of ray between origin and intersection


##############################
### Vector maths functions ###
##############################

cdef inline vector_t transform_c(transform_t t, vector_t p):
    cdef vector_t out
    out.x = p.x*t.m00 + p.y*t.m01 + p.z*t.m02 + t.tx
    out.y = p.x*t.m10 + p.y*t.m11 + p.z*t.m12 + t.ty
    out.z = p.x*t.m20 + p.y*t.m21 + p.z*t.m22 + t.tz
    return out

cdef inline vector_t set_v(vector_t v, object O):
    v.x = O[0]
    v.y = O[1]
    v.z = O[2]
    return v

cdef inline double sep(vector_t p1, vector_t p2):
    return sqrt((p2.x-p1.x)**2 + (p2.y-p1.y)**2 + (p2.z-p1.z)**2)

cdef inline vector_t multvv(vector_t a, vector_t b):
    cdef vector_t out
    out.x = a.x*b.x
    out.y = a.y*b.y
    out.z = a.z*b.z
    return out

cdef inline vector_t multvs(vector_t a, double b):
    cdef vector_t out
    out.x = a.x*b
    out.y = a.y*b
    out.z = a.z*b
    return out

cdef inline vector_t addvv(vector_t a, vector_t b):
    cdef vector_t out
    out.x = a.x+b.x
    out.y = a.y+b.y
    out.z = a.z+b.z
    return out


##################################
### Python extension types
##################################

cdef class Transform:
    cdef transform_t trans
    
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


cdef class Intersection:
    cdef intersection_t inter
    
    property face_idx:
        def __get__(self):
            return self.inter.face_idx
        
        def __set__(self, unsigned int i):
            self.inter.face_idx = i
            
    property point:
        def __get__(self):
            return (self.inter.point.x,
                    self.inter.point.y,
                    self.inter.point.z)
        def __set__(self, v):
            self.inter.point.x = v[0]
            self.inter.point.y = v[1]
            self.inter.point.z = v[2]
            
    property dist:
        def __get__(self):
            return self.inter.dist
        def __set__(self, double d):
            self.inter.dist = d


cdef class Ray:
    cdef ray_t ray
    
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
            return (self.ray.normals.x,self.ray.normals.y,self.ray.normals.z)
        
        def __set__(self, v):
            self.ray.normals.x = v[0]
            self.ray.normals.y = v[1]
            self.ray.normals.z = v[2]
            
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
    cdef ray_t *rays
    cdef readonly unsigned long n_rays, max_size
    cdef public RayCollection parent
    
    def __cinit__(self, size_t max_size):
        self.rays = <ray_t*>malloc(max_size*sizeof(ray_t))
        self.n_rays = 0
        self.max_size = max_size
        
    def __dealloc__(self):
        free(self.rays)
        
    cdef add_ray_c(self, ray_t r):
        self.rays[self.n_rays] = r
        self.n_rays += 1
        
    def add_ray(self, Ray r):
        if self.n_rays == self.max_size-1:
            raise ValueError("RayCollection is full up")
        self.rays[self.n_rays] = r.ray
        self.n_rays += 1
        
    def add_ray_list(self, list rays):
        cdef int i
        if self.n_rays + len(rays) >= self.max_size:
            raise IndexError("RayCollection size is too small to hold this many rays")
        for i in xrange(len(rays)):
            if not isinstance(rays[i], Ray):
                raise TypeError("ray list contains non-Ray instance at index %d"%i)
        for i in xrange(len(rays)):
            self.add_ray_c((<Ray>rays[i]).ray)
        
    def clear_ray_list(self):
        self.n_rays = 0
        
    def get_ray_list(self):
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
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        self.rays[idx] = r.ray
    
    
cdef class Face(object):
    cdef public object owner
    cdef public char *name
    cdef public double tolerance
    cdef public int idx #index in the global face list
    
    def __cinit__(self):
        self.name = "base Face class"
        self.tolerance = 0.0001
    
    cdef vector_t intersect_c(self, vector_t p1, vector_t p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        return p1
    
    def intersect(self, p1, p2):
        cdef:
            vector_t p1_, p2_, p_i
        
        set_v(p1_, p1)
        set_v(p2_, p2)
        p_i = self.intersect_c(p1_, p2_)
        return (p_i.x, p_i.y, p_i.z)

    cdef vector_t compute_normal_c(self, vector_t p):
        return p
    
    def compute_normal(self, p):
        cdef vector_t p_, n
        n = self.compute_normal_c(p_)
        return (n.x, n.y, n.z)
    
    cdef ray_t eval_child_ray_c(self, ray_t old_ray, int ray_idx, 
                                vector_t p):
        return old_ray
    
    def eval_child_ray(self, Ray old_ray, ray_idx, point):
        cdef:
            vector_t p
            Ray out=Ray()
            unsigned int idx
        
        set_v(p, point)
        out.ray = self.eval_child_ray_c(old_ray.ray, ray_idx, p)
        return out
        


cdef class FaceList(object):
    """A group of faces which share a transform"""
    cdef transform_t trans
    cdef transform_t inv_trans
    cdef public list faces
    
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
        
     
    cdef intersection_t intersect_c(self, vector_t P1, vector_t P2, double max_length):
        """Finds the face with the nearest intersection
        point, for the ray defined by the two input points,
        P1 and P2 (in global coords).
        """
        cdef vector_t p1, p2, point
        cdef list faces
        cdef double d, dist=INFINITY
        cdef Face nearest=None
        cdef intersection_t inter
        
        p1 = transform_c(self.trans, P1)
        p2 = transform_c(self.trans, P2)
        
        faces = self.faces
        
        for i in xrange(len(faces)):
            point = (<Face>(faces[i])).intersect_c(p1, p2)
            d = sep(p1, point)
            if 0.0 < d < dist:
                dist = d
                nearest = <Face>(faces[i])
        
        inter.face_idx = nearest.idx
        inter.point = transform_c(self.inv_trans, point)
        inter.dist = dist/sep(p1,p2)
        return inter
    
    def intersect(self, P1, P2, double max_length):
        cdef vector_t P1_, P2_
        cdef intersection_t inter
        cdef Intersection i2=Intersection()
        
        set_v(P1_, P1)
        set_v(P2_, P2)
        inter = self.intersect_c(P1_, P2_, max_length)
        i2.inter = inter
        return i2

##################################
### Python module functions
##################################

cdef RayCollection trace_segment_c(RayCollection rays, 
                                    list face_sets, 
                                    list all_faces):
    cdef:
        FaceList face_set #a FaceList
        int size, i
        vector_t P1, P2, normal
        float max_length
        intersection_t inter, inter2
        int n_sets=len(face_sets)
        ray_t ray
   
    #need to allocate the output rays here 
    new_rays = RayCollection(rays.n_rays)
   
    for i in range(size):
        ray = rays.rays[i]
        P1 = ray.origin
        P2 = addvv(P1, multvs(ray.direction, max_length))
        inter = (<FaceList>(face_sets[0])).intersect_c(P1, P2, max_length)
        dist = sep(inter.point, P1)
        for j in xrange(n_sets-1):
            face_set = face_sets[j+1]
            inter2 = (<FaceList>face_set).intersect_c(P1, P2, max_length)
            dist2 = sep(inter2.point, P1)
            if dist2 < dist:
                dist = dist2
                inter = inter2
                
        if dist <= INFINITY:
            face = all_faces[inter.face_idx]
            #evaluate new ray
            new_rays.add_ray_c(face.eval_child_ray_c(ray, i, 
                                                    inter.point))
    return new_rays


def transform(Transform t, p):
    cdef vector_t p1, p2
    assert isinstance(t, Transform)
    assert len(p)==3
    p1.x = p[0]
    p1.y = p[1]
    p1.z = p[2]
    p2 = transform_c(t.trans, p1)
    return (p2.x, p2.y, p2.z)
    


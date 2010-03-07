#!/bin/env python

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

cdef extern from "math.h":
    double sqrt(double arg)
    double fabs(double arg)
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

cdef inline vector_t rotate_c(transform_t t, vector_t p):
    cdef vector_t out
    out.x = p.x*t.m00 + p.y*t.m01 + p.z*t.m02
    out.y = p.x*t.m10 + p.y*t.m11 + p.z*t.m12
    out.z = p.x*t.m20 + p.y*t.m21 + p.z*t.m22
    return out

cdef inline vector_t set_v(vector_t v, object O):
    v.x = O[0]
    v.y = O[1]
    v.z = O[2]
    return v

def py_set_v(O):
    cdef vector_t v_
    print "before", (v_.x, v_.y, v_.z)
    v_ = set_v(v_, O)
    print "after", (v_.x, v_.y, v_.z)
    return (v_.x, v_.y, v_.z)

cdef inline double sep_(vector_t p1, vector_t p2):
    cdef double a,b
    a = (p2.x-p1.x)
    b = (p2.y-p1.y)
    c = (p2.z-p1.z)
    return sqrt((a*a) + (b*b) + (c*c))

def sep(a, b):
    cdef vector_t a_, b_
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
    return sep_(a_, b_)

cdef inline vector_t invert_(vector_t v):
    v.x = -v.x
    v.y = -v.y
    v.z = -v.z
    return v

def invert(v):
    cdef vector_t v_
    v_ = set_v(v_, v)
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
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
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
    a_ = set_v(a_, a)
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
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
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
    a_ = set_v(a_, a)
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
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
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
    a_ = set_v(a_, a)
    c_ = subvs_(a_, b)
    return (c_.x, c_.y, c_.z)

cdef inline double mag_(vector_t a):
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z)

def mag(a):
    cdef vector_t a_
    a_ = set_v(a_, a)
    return mag_(a_)

cdef inline double mag_sq_(vector_t a):
    return a.x*a.x + a.y*a.y + a.z*a.z

def mag_sq(a):
    cdef vector_t a_
    a_ = set_v(a_, a)
    return mag_sq_(a_)

cdef inline double dotprod_(vector_t a, vector_t b):
    return a.x*b.x + a.y*b.y + a.z*b.z

def dotprod(a, b):
    cdef vector_t a_, b_
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
    return dotprod_(a_,b_)

cdef inline vector_t cross_(vector_t a, vector_t b):
    cdef vector_t c
    c.x = a.y*b.z - a.z*b.y
    c.y = a.z*b.x - a.x*b.z
    c.z = a.x*b.y - a.y*b.x
    return c

def cross(a, b):
    cdef vector_t a_, b_, c_
    a_ = set_v(a_, a)
    b_ = set_v(b_, b)
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
    a_ = set_v(a_, a)
    a_ = norm_(a_)
    return (a_.x, a_.y, a_.z)

cdef ray_t convert_to_sp(ray_t ray, vector_t normal):
    """Project the E-field components of a given ray
    onto the S- and P-polarisations defined by the 
    surface normal
    """
    cdef:
        vector_t E2_vector, E1_vector, v, S_vector, P_vector
        complex_t S_amp, P_amp, E1_amp, E2_amp
        double A, B
    
    E1_amp = ray.E1_amp
    E2_amp = ray.E2_amp
    E1_vector = ray.E_vector
    E2_vector = cross_(ray.direction, E1_vector)
    v = cross_(ray.direction, normal)
    if v.x==0. and v.y==0. and v.z==0:
        S_vector = norm_(E1_vector)
    else:
        S_vector = norm_(v)
    v = cross_(ray.direction, S_vector)
    P_vector = norm_(v)
    
    A = dotprod_(E1_vector,S_vector)
    B = dotprod_(E2_vector, S_vector)
    
    S_amp.real = E1_amp.real*A + E2_amp.real*B
    S_amp.imag = E1_amp.imag*A + E2_amp.imag*B
    
    A = dotprod_(E1_vector, P_vector)
    B = dotprod_(E2_vector, P_vector)
    
    P_amp.real = E1_amp.real*A + E2_amp.real*B
    P_amp.imag = E1_amp.imag*A + E2_amp.imag*B
    
    ray.E_vector = S_vector
    ray.E1_amp = S_amp
    ray.E2_amp = P_amp
    return ray

def Convert_to_SP(Ray ray, normal):
    cdef vector_t n
    cdef ray_t r
    n = set_v(n, normal)
    r = convert_to_sp(ray.ray, n)
    out = Ray()
    out.ray = r
    return out

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


cdef class Intersection:
    
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
        
    cdef add_ray_c(self, ray_t r):
        self.rays[self.n_rays] = r
        self.n_rays += 1
        
    def reset_length(self):
        cdef int i
        for i in xrange(self.n_rays):
            self.rays[i].length = INFINITY
        
    def add_ray(self, Ray r):
        if self.n_rays == self.max_size:
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
    
    
cdef class InterfaceMaterial(object):
    """Abstract base class for objects describing
    the materials characterics of a Face
    """
    
    cdef ray_t eval_child_ray_c(self, ray_t old_ray, 
                                unsigned int ray_idx, 
                                vector_t p, vector_t normal):
        return old_ray
    
    def eval_child_ray(self, Ray old_ray, ray_idx, point, 
                        normal):
        cdef:
            vector_t p, n
            Ray out=Ray()
            unsigned int idx
        
        p = set_v(p, point)
        n = set_v(n, normal)
        out.ray = self.eval_child_ray_c(old_ray.ray, ray_idx, 
                                        p, n)
        return out
    
    
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
            self.material = PECMaterial()
        self.invert_normal = int(kwds.get('invert_normal', 0))
        
    
    cdef intersection_t intersect_c(self, vector_t p1, vector_t p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        cdef intersection_t inter
        return inter
    
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
            intersection_t inter
            Intersection py_inter
        
        p1_ = set_v(p1_, p1)
        p2_ = set_v(p2_, p2)
        inter = self.intersect_c(p1_, p2_)
        py_inter = Intersection()
        py_inter.inter = inter
        return py_inter

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
        
     
    cdef intersection_t intersect_c(self, vector_t P1, vector_t P2, double max_length):
        """Finds the face with the nearest intersection
        point, for the ray defined by the two input points,
        P1 and P2 (in global coords).
        """
        cdef vector_t p1, p2, point
        cdef list faces
        cdef unsigned int i
        cdef intersection_t inter, nearest
        
        p1 = transform_c(self.inv_trans, P1)
        p2 = transform_c(self.inv_trans, P2)
        #print "LOCAL", p1, p2
        faces = self.faces
        
        nearest.dist = INFINITY
        
        for i in xrange(len(faces)):
            inter = (<Face>(faces[i])).intersect_c(p1, p2)
            if 0.0 < inter.dist < nearest.dist:
                nearest = inter
        
        nearest.point = transform_c(self.trans, nearest.point)
        
        #print "INTER", nearest.point, nearest.dist, nearest.face_idx
        return nearest
    
    def intersect(self, P1, P2, double max_length):
        cdef vector_t P1_, P2_
        cdef intersection_t inter
        cdef Intersection i2=Intersection()
        
        P1_ = set_v(P1_, P1)
        P2_ = set_v(P2_, P2)
        inter = self.intersect_c(P1_, P2_, max_length)
        i2.inter = inter
        return i2
    
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
        p = set_v(p, point)
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
        unsigned int size, i, j, nearest_set=0
        vector_t P1, P2, normal, point
        intersection_t inter, nearest
        int n_sets=len(face_sets)
        ray_t ray, new_ray
        RayCollection new_rays
   
    #need to allocate the output rays here 
    new_rays = RayCollection(rays.n_rays)
    
    for i in range(rays.n_rays):
        ray = rays.rays[i]
        P1 = ray.origin
        P2 = addvv_(P1, multvs_(ray.direction, max_length))
        #print "points", P1, P2
        nearest = (<FaceList>(face_sets[0])).intersect_c(P1, P2, max_length)
        for j in xrange(n_sets-1):
            face_set = face_sets[j+1]
            inter = (<FaceList>face_set).intersect_c(P1, P2, max_length)
            if inter.dist < nearest.dist:
                nearest = inter
                nearest_set = j+1
                
        if nearest.dist < INFINITY:
            #print "GET FACE", nearest.face_idx, len(all_faces)
            face = all_faces[nearest.face_idx]
            #evaluate new ray
            ray.end_face_idx = nearest.face_idx
            rays.rays[i].length = nearest.dist
            #print "ray length", ray.length
            point = nearest.point
            normal = (<FaceList>(face_sets[nearest_set])).compute_normal_c(face, point)
            #print "s normal", normal
            new_ray = (<InterfaceMaterial>(face.material)).eval_child_ray_c(ray, i, 
                                                    point,
                                                    normal
                                                    )
            new_ray.length = max_length
            new_rays.add_ray_c(new_ray)
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
    
    
    
cdef class PECMaterial(InterfaceMaterial):
    """Simulates a Perfect Electrical Conductor
    """
    cdef ray_t eval_child_ray_c(self,
                                ray_t in_ray, 
                                unsigned int idx, 
                                vector_t point,
                                vector_t normal):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected
            ray_t sp_ray
            complex_t cpx
            double cosTheta
        
        normal = norm_(normal)
        sp_ray = convert_to_sp(in_ray, normal)
        cosTheta = dotprod_(normal, in_ray.direction)
        cosThetaNormal = multvs_(normal, cosTheta)
        reflected = subvv_(in_ray.direction, multvs_(cosThetaNormal, 2))
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.direction = reflected
        sp_ray.E1_amp.real = -sp_ray.E1_amp.real
        sp_ray.E1_amp.imag = -sp_ray.E1_amp.imag
        sp_ray.E2_amp.real = -sp_ray.E2_amp.real
        sp_ray.E2_amp.imag = -sp_ray.E2_amp.imag
        sp_ray.parent_idx = idx
        return sp_ray
    
    
cdef class DielectricMaterial(InterfaceMaterial):
    """Simulates Fresnel reflection and refraction at a
    normal dielectric interface
    """
    
    def __cinit__(self, **kwds):
        self.n_inside = kwds.get('n_inside', 1.5)
        self.n_outside = kwds.get('n_outside', 1.0)
    
    property n_inside:
        def __get__(self):
            return complex(self.n_inside_.real, self.n_inside_.imag)
        
        def __set__(self, v):
            v = complex(v)
            self.n_inside_.real = v.real
            self.n_inside_.imag = v.imag
            
    property n_outside:
        def __get__(self):
            return complex(self.n_outside_.real, self.n_outside_.imag)
        
        def __set__(self, v):
            v = complex(v)
            self.n_outside_.real = v.real
            self.n_outside_.imag = v.imag
            
    cdef ray_t eval_child_ray_c(self,
                                ray_t in_ray, 
                                unsigned int idx, 
                                vector_t point,
                                vector_t normal):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2
            ray_t sp_ray
            complex_t cpx
            double cosTheta, n1, n2, N2, cos1
            double N2cosTheta, N2_sin2, tan_mag_sq, c2
            double cos2, Two_n1_cos1, aspect, T_p, T_s
            int flip
            
        normal = norm_(normal)
        in_ray.direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray, normal)
        cosTheta = dotprod_(normal, in_ray.direction)
        cos1 = fabs(cosTheta)
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside_.real
            n2 = self.n_inside_.real
            sp_ray.refractive_index = self.n_inside_
            flip = 1
            print "out to in", n1, n2
        else:
            n1 = self.n_inside_.real
            n2 = self.n_outside_.real
            sp_ray.refractive_index = self.n_outside_
            flip = -1
            print "in to out", n1, n2
            
        N2 = (n2/n1)**2
        N2cosTheta = N2*cos1
        
        N2_sin2 = (cosTheta*cosTheta) + (N2 - 1)
        print "TIR", N2_sin2, cosTheta, N2, cos1
        print (normal.x, normal.y, normal.z), in_ray.direction
        cosThetaNormal = multvs_(normal, cosTheta)
        if N2_sin2 < 0.0:
            #total internal reflection
            reflected = subvv_(in_ray.direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INFINITY
            sp_ray.E1_amp.real *= -1
            sp_ray.E1_amp.imag *= -1
            sp_ray.E2_amp.real *= -1
            sp_ray.E2_amp.imag *= -1
            sp_ray.parent_idx = idx
            return sp_ray
        else:
            #normal transmission            
            tangent = subvv_(in_ray.direction, cosThetaNormal)
            tg2 = multvs_(tangent, n1/n2)
            tan_mag_sq = mag_sq_(tg2)
            c2 = sqrt(1-tan_mag_sq)
            transmitted = subvv_(tg2, multvs_(normal, c2*flip))
            
            cos2 = fabs(dotprod_(transmitted, normal))
            Two_n1_cos1 = (2*n1)*cos1
            aspect = sqrt(cos2/cos1) * Two_n1_cos1
            
            #Fresnel equations for transmission
            T_p = aspect / ( n2*cos1 + n1*cos2 )
            T_s = aspect / ( n2*cos2 + n1*cos1 )
            #print "T_s", T_s, "T_p", T_p
            
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = transmitted
            sp_ray.length = INFINITY
            sp_ray.E1_amp.real *= T_s
            sp_ray.E1_amp.imag *= T_s
            sp_ray.E2_amp.real *= T_p
            sp_ray.E2_amp.imag *= T_p
            sp_ray.parent_idx = idx
            return sp_ray
            return sp_ray
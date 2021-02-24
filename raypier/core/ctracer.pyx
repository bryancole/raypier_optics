#!/bin/env python

#cython: boundscheck=False
#cython: nonecheck=False
#cython: cdivision=True

cdef extern from "math.h":
    double M_PI
    double sqrt(double arg) nogil
    double fabs(double arg) nogil
    double atan2(double y, double x) nogil
    double atan(double arg) nogil
    double sin(double arg) nogil
    double cos(double arg) nogil
    #double INFINITY
    
cdef extern from "float.h":
    double DBL_MAX

cdef:
    INF=(DBL_MAX+DBL_MAX)
    double root2 = sqrt(2.0)
    public unsigned int REFL_RAY=1<<0
    public unsigned int GAUSSLET=1<<1
    public unsigned int PARABASAL=1<<2
    
from libc.stdlib cimport malloc, free, realloc

cdef extern from "stdlib.h" nogil:
    void *memcpy(void *str1, void *str2, size_t n)

import time
import numpy as np
cimport numpy as np_

cdef:
    int NPARA = 6
    

ray_dtype = np.dtype([('origin', np.double, (3,)),
                        ('direction', np.double, (3,)),
                        ('normal', np.double, (3,)),
                        ('E_vector', np.double, (3,)),
                        ('refractive_index', np.complex128),
                        ('E1_amp', np.complex128),
                        ('E2_amp', np.complex128),
                        ('length', np.double),
                        ('phase', np.double),
                        ('accumulated_path', np.double),
                        ('wavelength_idx', np.uint32),
                        ('parent_idx', np.uint32),
                        ('end_face_idx', np.uint32),
                        ('ray_type_id', np.uint32) #No point using a smaller type here as it'll probably get padded by the compiler
                        ])
                        
    
GAUSSLET_ = GAUSSLET
PARABASAL_ = PARABASAL

para_dtype = np.dtype([('origin', np.double, (3,)),
                        ('direction', np.double, (3,)),
                        ('normal', np.double, (3,)),
                        ('length', np.double),
                        ])


gausslet_dtype = np.dtype([
                    ('base_ray', ray_dtype),
                    ('para_rays', para_dtype, (NPARA,))
                    ])


##############################
### Vector maths functions ###
##############################

cdef inline vector_t transform_c(transform_t t, vector_t p) nogil:
    cdef vector_t out
    out.x = p.x*t.m00 + p.y*t.m01 + p.z*t.m02 + t.tx
    out.y = p.x*t.m10 + p.y*t.m11 + p.z*t.m12 + t.ty
    out.z = p.x*t.m20 + p.y*t.m21 + p.z*t.m22 + t.tz
    return out

cdef inline vector_t rotate_c(transform_t t, vector_t p) nogil:
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

cdef inline double sep_(vector_t p1, vector_t p2) nogil:
    cdef double a,b
    a = (p2.x-p1.x)
    b = (p2.y-p1.y)
    c = (p2.z-p1.z)
    return sqrt((a*a) + (b*b) + (c*c))

def sep(a, b):
    cdef vector_t a_ = set_v(a), b_ = set_v(b)
    return sep_(a_, b_)

cdef inline vector_t invert_(vector_t v) nogil:
    v.x = -v.x
    v.y = -v.y
    v.z = -v.z
    return v

def invert(v):
    cdef vector_t v_ = set_v(v)
    v_ = invert_(v_)
    return (v_.x, v_.y, v_.z)

cdef inline vector_t multvv_(vector_t a, vector_t b) nogil:
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

cdef inline vector_t multvs_(vector_t a, double b) nogil:
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

cdef inline vector_t addvv_(vector_t a, vector_t b) nogil:
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
    
cdef inline vector_t addvs_(vector_t a, double b) nogil:
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

cdef inline vector_t subvv_(vector_t a, vector_t b) nogil:
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

cdef inline vector_t subvs_(vector_t a, double b) nogil:
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

cdef inline double mag_(vector_t a) nogil:
    return sqrt(a.x*a.x + a.y*a.y + a.z*a.z)

def mag(a):
    cdef vector_t a_
    a_ = set_v(a)
    return mag_(a_)

cdef inline double mag_sq_(vector_t a) nogil:
    return a.x*a.x + a.y*a.y + a.z*a.z

def mag_sq(a):
    cdef vector_t a_
    a_ = set_v(a)
    return mag_sq_(a_)

cdef inline double dotprod_(vector_t a, vector_t b) nogil:
    return a.x*b.x + a.y*b.y + a.z*b.z

def dotprod(a, b):
    cdef vector_t a_, b_
    a_ = set_v(a)
    b_ = set_v(b)
    return dotprod_(a_,b_)

cdef inline vector_t cross_(vector_t a, vector_t b) nogil:
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

cdef vector_t norm_(vector_t a) nogil:
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
    
    
cdef class GaussletCollectionIterator:        
    def __cinit__(self, GaussletCollection rays):
        self.rays = rays
        self.counter = 0
        
    def __iter__(self):
        return self
    
    def __next__(self):
        cdef Gausslet ray=Gausslet.__new__(Gausslet)
        if self.counter >= self.rays.n_rays:
            raise StopIteration
        ray.gausslet = self.rays.rays[self.counter]
        self.counter += 1
        return ray


cdef class ParabasalRay:
    def __cinit__(self, **kwds):
        self.max_length = 1000.0
        for k in kwds:
            setattr(self, k, kwds[k])
            
    def __repr__(self):
        return "Parabasal Ray(o=%s, d=%s)"%(str(self.origin),
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
            
    property normal:
        """normal vector for the face which created this ray"""
        def __get__(self):
            return (self.ray.normal.x,self.ray.normal.y,self.ray.normal.z)
        
        def __set__(self, v):
            self.ray.normal.x = v[0]
            self.ray.normal.y = v[1]
            self.ray.normal.z = v[2]
            
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
            if self.ray.length > self.max_length:
                length = self.max_length
            else:
                length = self.ray.length
            end = addvv_(self.ray.origin, multvs_(self.ray.direction, 
                                    length))
            return (end.x, end.y, end.z)


cdef class Ray:
    """ Ray - a wrapper around the ray_t C-structure.
    
    The Ray extension class exists mainly as a convenience for manipulation of single or small numbers of rays 
    from python. Large numbers of rays are more efficiently handled as either RayCollection objects, created in the
    tracing process, or as numpy arrays with the 'ray_dtype' dtype.
    """
    
    def __cinit__(self, **kwds):
        self.max_length = 1000.0
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
            
    property normal:
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
            
    property phase:
        """An additional phase-factor for the ray. At present, this handles the 'grating phase' factor
        generated by diffraction gratings. All other material surfaces leave this unchanged"""
        def __get__(self):
            return self.ray.phase
        
        def __set__(self, double v):
            self.ray.phase = v
            
    property accumulated_path:
        """The total *optical* path up to the start-point of this ray."""
        def __get__(self):
            return self.ray.accumulated_path
        
        def __set__(self, double v):
            self.ray.accumulated_path = v
            
    property wavelength_idx:
        """The wavelength of the ray in vacuum, in microns"""
        def __get__(self):
            return self.ray.wavelength_idx
        
        def __set__(self, unsigned int v):
            self.ray.wavelength_idx = v
            
    property termination:
        """the end-point of the ray (read only)
        """
        def __get__(self):
            cdef vector_t end
            cdef float length
            if self.ray.length > self.max_length:
                length = self.max_length
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
        
        def __set__(self, unsigned int v):
            self.ray.parent_idx = v
            
    property end_face_idx:
        """Index of the terminating face, in the global
        face list (created for each tracing operation)
        """
        def __get__(self):
            return self.ray.end_face_idx
        
        def __set__(self, unsigned int v):
            self.ray.end_face_idx = v
            
    property ray_type_id:
        """Used to distinguish rays created by reflection vs transmission or some other mechanism.
        Transmission->0, Reflection->1"""
        def __get__(self):
            return self.ray.ray_type_id
        
        def __set__(self, unsigned int v):
            self.ray.ray_type_id = v
            
    property power:
        """Optical power for the ray"""
        def __get__(self):
            return ray_power_(self.ray)
        
    property amplitude:
        """E field amplitude"""
        def __get__(self):
            cdef ray_t ray=self.ray
            return sqrt(ray.E1_amp.real**2 + ray.E1_amp.imag**2 +\
                        ray.E2_amp.real**2 + ray.E2_amp.imag**2)
            
    property jones_vector:
        """Jones polarisation vector expressed as a tuple (alpha, beta)
        where alpha and beta are complex"""
        def __get__(self):
            cdef double amp = self.amplitude
            return (self.E1_amp/amp, self.E2_amp/amp)
        
    property E_left:
        def __get__(self):
            E1_amp = self.E1_amp
            E2_amp = self.E2_amp
            return (E1_amp + 1.0j*E2_amp)/sqrt(2)
        
    property E_right:
        def __get__(self):
            E1_amp = self.E1_amp
            E2_amp = self.E2_amp
            return (E1_amp - 1.0j*E2_amp)/sqrt(2)
        
    property ellipticity:
        """Provide the ratio of power in the RH circular polarisation
        to the LH circular polarisation. A value of zero indicates 
        linear polarisation. +1 indicate RH polarisation, -1 is
        LH polarisation. Or maybe the other way round."""
        def __get__(self):
            alpha, beta = self.jones_vector
            R = (alpha - 1j*beta)/root2
            L = (alpha + 1j*beta)/root2
            PR = R.real**2 + R.imag**2
            PL = L.real**2 + L.imag**2
            return (PR - PL)/(PR + PL)
        
    property major_minor_axes:
        """Find the vector giving the major and minor axes of the polarisation ellipse.
        For fully circularly polarised light, the current E_vector will be
        returned"""
        def __get__(self):
            cdef:
                vector_t E1_vector, E2_vector
                
            E2_vector = norm_(cross_(self.ray.direction, self.ray.E_vector))
            E1_vector = norm_(cross_(E2_vector, self.ray.direction))
            
            E1 = sqrt(self.ray.E1_amp.real**2 + self.ray.E1_amp.imag**2)
            E2 = sqrt(self.ray.E2_amp.real**2 + self.ray.E2_amp.imag**2)
            phi1 = atan2(self.ray.E1_amp.imag, self.ray.E1_amp.real)
            phi2 = atan2(self.ray.E2_amp.imag, self.ray.E2_amp.real)
            
            A = E1*E1*sin(2*phi1) + E2*E2*sin(2*phi2)
            B = E1*E1*cos(2*phi1) + E2*E2*cos(2*phi2)
            
            mag = sqrt(A*A + B*B)
            t1 = atan( (B - mag)/A )
            t2 = atan( (B + mag)/A )
            
            axis1 = addvv_(multvs_(E1_vector, E1*cos(phi1+t1)),
                           multvs_(E2_vector, E2*cos(phi2+t1)))
            
            axis2 = addvv_(multvs_(E1_vector, E1*cos(phi1+t2)),
                           multvs_(E2_vector, E2*cos(phi2+t2)))
            
            return ( (axis1.x, axis1.y, axis1.z), 
                     (axis2.x, axis2.y, axis2.z) )
            
        
    def project_E(self, *axis):
        """Rotate the E_vector onto the given axis, projecting
        E1_amp and E2_amp as necessary."""
        cdef:
            complex_t E1, E2
            vector_t v=set_v(axis)
            vector_t E_vector, E_vector2bar, E_vector2, direction
            double A, B
            
        E1 = self.ray.E1_amp
        E2 = self.ray.E2_amp
        direction = self.ray.direction
        
        #initial E_vector
        E_vector = norm_(cross_(cross_(direction, self.ray.E_vector), direction))
        #new E_vector
        E_vector2bar = norm_(cross_(direction, v))
        E_vector2 = norm_(cross_(E_vector2bar, direction))
        
        A = dotprod_(E_vector, E_vector2)
        B = dotprod_(E_vector, E_vector2bar)
        
        #print "A,B,AB^2:", A, B, A**2 + B**2
        
        self.ray.E1_amp.real = -(E2.real*B - E1.real*A)
        self.ray.E2_amp.real = (E2.real*A + E1.real*B)
        
        self.ray.E1_amp.imag = -(E2.imag*B - E1.imag*A)
        self.ray.E2_amp.imag = E2.imag*A + E1.imag*B
        
        self.ray.E_vector = E_vector2
        
        
cdef class Gausslet:
    def __cinit__(self, **kwds):
        for k in kwds:
            setattr(self, k, kwds[k])
            
    def __repr__(self):
        return "Gausslet(o=%s, d=%s)"%(str(self.base_ray.origin),
                                            str(self.base_ray.direction))
        
    property base_ray:
        def __get__(self):
            cdef:
                Ray out = Ray()
                
            out.ray = self.gausslet.base_ray
            return out
        
        def __set__(self, Ray r):
            self.gausslet.base_ray = r.ray
            
    property parabasal_rays:
        def __get__(self):
            cdef:
                ParabasalRay p
                int i
                tuple out
                
            out = tuple(ParabasalRay() for i in range(6))
            for i in range(6):
                p = out[i]
                p.ray = self.gausslet.para[i]
            return out
        
        def __set__(self, tuple paras):
            for i in range(6):
                self.gausslet.para[i] = paras[i].ray
                
                
cdef class RayArrayView:
    """An abstract class to provide the API for ray_t member access from python / numpy"""
    cdef void set_ray_c(self, unsigned long i, ray_t ray):
        return
    
    cdef ray_t get_ray_c(self, unsigned long i):
        cdef:
            ray_t r
        return r
    
    cdef unsigned long get_n_rays(self):
        return 0
    
    def __getitem__(self, size_t idx):
        cdef Ray r
        if idx >= self.get_n_rays():
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        r = Ray()
        r.ray = self.get_ray_c(idx)
        return r
    
    def __setitem__(self, size_t idx, Ray r):
        if idx >= self.get_n_rays():
            raise IndexError("Attempting to set index %d from a size %d array"%(idx, self.n_rays))
        self.set_ray_c(idx, r.ray)
        
    def get_ray_list(self):
        """Returns the contents of this RayCollection as a list of Rays
        """
        cdef size_t i
        cdef list ray_list = []
        cdef Ray r
        for i in range(self.get_n_rays()):
            r = Ray()
            r.ray = self.get_ray_c(i)
            ray_list.append(r)
        return ray_list
    
    property origin:
        def __get__(self):
            cdef:
                unsigned long n=self.get_n_rays()
                np_.ndarray out = np.empty((n,3), dtype='d')
                size_t i
                vector_t v
            for i in xrange(n):
                v = self.get_ray_c(i).origin
                out[i,0] = v.x
                out[i,1] = v.y
                out[i,2] = v.z
            return out
        
    property direction:
        def __get__(self):
            cdef:
                unsigned long n=self.get_n_rays()
                np_.ndarray out = np.empty((n,3), dtype='d')
                size_t i
                vector_t v
            for i in xrange(n):
                v = self.get_ray_c(i).direction
                out[i,0] = v.x
                out[i,1] = v.y
                out[i,2] = v.z
            return out
        
    property normal:
        def __get__(self):
            cdef:
                unsigned long n=self.get_n_rays()
                np_.ndarray out = np.empty((n,3), dtype='d')
                size_t i
                vector_t v
            for i in xrange(n):
                v = self.get_ray_c(i).normal
                out[i,0] = v.x
                out[i,1] = v.y
                out[i,2] = v.z
            return out
        
    property E_vector:
        def __get__(self):
            cdef:
                unsigned long n=self.get_n_rays()
                np_.ndarray out = np.empty((n,3), dtype='d')
                size_t i
                vector_t v
            for i in xrange(n):
                v = self.get_ray_c(i).E_vector
                out[i,0] = v.x
                out[i,1] = v.y
                out[i,2] = v.z
            return out
        
    property refractive_index:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.complex128)
                np_.ndarray rr,ii
                size_t i
                complex_t n
            rr = out.real
            ii = out.imag
            for i in xrange(N):
                n = self.get_ray_c(i).refractive_index
                rr[i] = n.real
                ii[i] = n.imag
            return out
        
    property E1_amp:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.complex128)
                np_.ndarray rr,ii
                size_t i
                complex_t n
            rr = out.real
            ii = out.imag
            for i in xrange(N):
                n = self.get_ray_c(i).E1_amp
                rr[i] = n.real
                ii[i] = n.imag
            return out
        
    property E2_amp:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.complex128)
                np_.ndarray rr,ii
                size_t i
                complex_t n
            rr = out.real
            ii = out.imag
            for i in xrange(N):
                n = self.get_ray_c(i).E2_amp
                rr[i] = n.real
                ii[i] = n.imag
            return out
        
    property length:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype='d')
                size_t i
                double v
            for i in xrange(N):
                v = self.get_ray_c(i).length
                out[i] = v
            return out
        
    property phase:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype='d')
                size_t i
                double v
            for i in xrange(N):
                v = self.get_ray_c(i).phase
                out[i] = v
            return out
        
    property accumulated_path:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty((N,), dtype='d' )
                size_t i
            for i in xrange(N):
                out[i] = self.get_ray_c(i).accumulated_path
            return out
        
    property wavelength_idx:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.uint32)
                size_t i
                double v
            for i in xrange(N):
                v = self.get_ray_c(i).wavelength_idx
                out[i] = v
            return out
        
    property parent_idx:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.uint32)
                size_t i
                unsigned int v
            for i in xrange(N):
                v = self.get_ray_c(i).parent_idx
                out[i] = v
            return out
        
    property end_face_idx:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.uint32)
                size_t i
                unsigned int v
            for i in xrange(N):
                v = self.get_ray_c(i).end_face_idx
                out[i] = v
            return out
        
    property ray_type_id:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty(N, dtype=np.uint32)
                size_t i
                unsigned int v
            for i in xrange(N):
                v = self.get_ray_c(i).ray_type_id
                out[i] = v
            return out
        
    property termination:
        def __get__(self):
            cdef:
                unsigned long N=self.get_n_rays()
                np_.ndarray out = np.empty((N,3), dtype='d')
                size_t i
                vector_t v
                ray_t r
            for i in xrange(N):
                r = self.get_ray_c(i)
                v = addvv_(r.origin, 
                           multvs_(r.direction, 
                                    r.length))
                out[i,0] = v.x
                out[i,1] = v.y
                out[i,2] = v.z
            return out
    
    

cdef class RayCollection(RayArrayView):
    """A list-like collection of ray_t objects.
    
    The RayCollection is the primary data-structure used in the ray-tracing operation. 
    
    The RayCollection is of variable length, in that it can grow as individual rays are added to it.
    Internally, the memory allocated to the array of ray_t structures is re-allocated to increase
    its capacity.
    """
    
    def __cinit__(self, size_t max_size):
        self.rays = <ray_t*>malloc(max_size*sizeof(ray_t))
        self.n_rays = 0
        self.max_size = max_size
        self._mtime = 0.0
        
    def __dealloc__(self):
        free(self.rays)
        
    def __len__(self):
        return self.n_rays
    
    cdef ray_t get_ray_c(self, unsigned long i):
        return self.rays[i]
    
    cdef void set_ray_c(self, unsigned long i, ray_t ray):
        self.rays[i] = ray
        
    cdef unsigned long get_n_rays(self):
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
        
    cdef void reset_length_c(self, double max_length):
        cdef:
            size_t i
        for i in range(self.n_rays):
            self.rays[i].length = max_length
            
    def reset_length(self, double max_length=INF):
        """Sets the length of all rays in this RayCollection to Infinity
        """
        self.reset_length_c(max_length)
        
    def add_ray(self, Ray r):
        """Adds the given Ray instance to this collection
        """
        self.add_ray_c(r.ray)
        
    def add_ray_list(self, list rays):
        """Adds the given list of Rays to this collection
        """
        cdef long int i
        for i in range(len(rays)):
            if not isinstance(rays[i], Ray):
                raise TypeError("ray list contains non-Ray instance at index %d"%i)
        for i in range(len(rays)):
            self.add_ray_c((<Ray>rays[i]).ray)
        
    def clear_ray_list(self):
        """Empties this RayCollection (by setting the count to zero)
        """
        self.n_rays = 0
        
    def __iter__(self):
        return RayCollectionIterator(self)
    
    def copy_as_array(self):
        """Returns the contents of this RayCollection as a numpy array
        (the data is always copied).
        """
        cdef np_.ndarray out = np.empty(self.n_rays, dtype=ray_dtype)
        memcpy(<np_.float64_t *>out.data, self.rays, self.n_rays*sizeof(ray_t))
        return out
    
    property wavelengths:
        def __get__(self):
            return np.asarray(self._wavelengths)
        
        def __set__(self, wl_list):
            self._wavelengths = np.ascontiguousarray(wl_list, dtype=np.double)
    
    cdef double get_mtime(self, unsigned long guard):
        cdef:
            double pmtime
            
        if self._parent is None:
            return self._mtime
        if guard == id(self):
            return self._mtime
        
        pmtime = self._parent.get_mtime(guard)
        if pmtime > self._mtime:
            self._neighbours = None
            return pmtime
        else:
            return self._mtime
        
    cdef void _eval_neighbours(self, int[:,:] pnb):
        cdef:
            unsigned long i
            int j, pidx, child_nb, rtype
            int[:,:] rmap
            unsigned int nparent=self._parent.n_rays
            
        if pnb is None:
            return
        
        rmap = np.full( (nparent, 2), -1 ,dtype=np.int32)
        nb = np.full( (self.n_rays, pnb.shape[1]), -1, dtype=np.int32)
        
        for i in range(self.n_rays):
            rtype = self.rays[i].ray_type_id & REFL_RAY
            pidx = self.rays[i].parent_idx
            rmap[pidx, rtype] = i
            
        for i in range(self.n_rays):
            rtype = self.rays[i].ray_type_id & REFL_RAY
            pidx = self.rays[i].parent_idx
            for j in range(pnb.shape[1]):
                child_nb = pnb[pidx,j]
                if child_nb >= 0:
                    nb[i,j] = rmap[child_nb,rtype]
                    
        self._neighbours = nb
        
        
    property neighbours:
        def __get__(self):
            if self._parent is None:
                if self._neighbours is None:
                    return None
                else:
                    return np.asarray(self._neighbours)
            else:
                pmtime = self._parent.get_mtime(id(self))
                if self._mtime >= pmtime:
                    if self._neighbours is not None:
                        return np.asarray(self._neighbours)
                self._eval_neighbours(self._parent.neighbours)
                self._mtime = time.monotonic()
                return np.asarray(self._neighbours)
                
        def __set__(self, int[:,:] nb):
            self._neighbours = nb
            self._mtime = time.monotonic()
            
    property parent:
        def __get__(self):
            return self._parent
        
        def __set__(self, RayCollection rc):
            self._parent = rc
            self._neighbours = None
            self._wavelengths = rc._wavelengths
    
    @classmethod
    def from_array(cls, np_.ndarray data):
        """Creates a new RayCollection from the given numpy array. The array
        dtype should be a ctracer.ray_dtype. The data is copied into the 
        RayCollection
        """
        cdef int size=data.shape[0]
        cdef RayCollection rc = RayCollection(size)
        assert data.dtype is ray_dtype
        data = np.ascontiguousarray(data)
        memcpy(rc.rays, <np_.float64_t *>data.data, size*sizeof(ray_t))
        rc.n_rays = size
        return rc
    
    
cdef class GaussletBaseRayView(RayArrayView):
    def __cinit__(self, GaussletCollection owner):
        self.owner = owner

    cdef void set_ray_c(self, unsigned long i, ray_t ray):
        self.owner.rays[i].base_ray = ray
    
    cdef ray_t get_ray_c(self, unsigned long i):
        return self.owner.rays[i].base_ray
    
    cdef unsigned long get_n_rays(self):
        return self.owner.n_rays
    
    def __len__(self):
        return self.owner.n_rays
    
    def copy_as_array(self):
        cdef:
            unsigned int i, N = self.get_n_rays()
            np_.ndarray out = np.empty((N,), dtype=ray_dtype)
            ray_t * rays ###A bug in Cython memory views means we cannot use ray_t[:]
            
        rays = <ray_t *>(out.data)
        for i in range(N):
            rays[i] = self.owner.rays[i].base_ray
        return out
    
    
cdef class GaussletCollection:
    """A list-like collection of ray_t objects.
    
    The RayCollection is the primary data-structure used in the ray-tracing operation. 
    
    The RayCollection is of variable length, in that it can grow as individual rays are added to it.
    Internally, the memory allocated to the array of ray_t structures is re-allocated to increase
    its capacity.
    """
    
    def __cinit__(self, size_t max_size):
        self.rays = <gausslet_t*>malloc(max_size*sizeof(gausslet_t))
        self.n_rays = 0
        self.max_size = max_size
        
    def __dealloc__(self):
        free(self.rays)
        
    def __len__(self):
        return self.n_rays
    
    property parent:
        def __get__(self):
            return self._parent
        
        def __set__(self, GaussletCollection gc):
            self._parent = gc
            self._wavelengths = gc._wavelengths
        
    cdef add_gausslet_c(self, gausslet_t r):
        if self.n_rays == self.max_size:
            if self.max_size == 0:
                self.max_size = 1
            else:
                self.max_size *= 2
            self.rays = <gausslet_t*>realloc(self.rays, self.max_size*sizeof(gausslet_t))
        self.rays[self.n_rays] = r
        self.n_rays += 1
        
    def add_gausslet(self, Gausslet r):
        """Adds the given Ray instance to this collection
        """
        self.add_gausslet_c(r.gausslet)
        
    def add_gausslet_list(self, list rays):
        """Adds the given list of Rays to this collection
        """
        cdef long int i
        for i in range(len(rays)):
            if not isinstance(rays[i], Gausslet):
                raise TypeError("ray list contains non-Gausslet instance at index %d"%i)
        for i in range(len(rays)):
            self.add_ray_c((<Gausslet>rays[i]).gausslet)
            
    cdef void reset_length_c(self, double max_length):
        cdef:
            size_t i, j
        for i in range(self.n_rays):
            self.rays[i].base_ray.length = max_length
            for j in range(6):
                self.rays[i].para[j].length = max_length
        
            
    def reset_length(self, double max_length=INF):
        """Sets the length of all rays in this RayCollection to Infinity
        """
        self.reset_length_c(max_length)
        
    def clear_ray_list(self):
        """Empties this RayCollection (by setting the count to zero)
        """
        self.n_rays = 0
        
    def get_gausslet_list(self):
        """Returns the contents of this RayCollection as a list of Rays
        """
        cdef size_t i
        cdef list ray_list = []
        cdef Gausslet r
        for i in range(self.n_rays):
            r = Gausslet()
            r.gausslet = self.rays[i]
            ray_list.append(r)
        return ray_list
    
    def __getitem__(self, size_t idx):
        cdef Gausslet r
        if idx >= self.n_rays:
            raise IndexError("Requested index %d from a size %d array"%(idx, self.n_rays))
        r = Gausslet()
        r.gausslet = self.rays[idx]
        return r
    
    def __setitem__(self, size_t idx, Gausslet r):
        if idx >= self.n_rays:
            raise IndexError("Attempting to set index %d from a size %d array"%(idx, self.n_rays))
        self.rays[idx] = r.gausslet
        
    def __iter__(self):
        return GaussletCollectionIterator(self)
    
    def copy_as_array(self):
        """Returns the contents of this RayCollection as a numpy array
        (the data is always copied).
        """
        cdef np_.ndarray out = np.empty(self.n_rays, dtype=gausslet_dtype)
        memcpy(<np_.float64_t *>out.data, self.rays, self.n_rays*sizeof(gausslet_t))
        return out
    
    def extend(self, GaussletCollection gc):
        self.extend_c(gc)
        
    cdef void extend_c(self, GaussletCollection gc):
        if (self.n_rays + gc.n_rays) > self.max_size:
            self.max_size = (self.n_rays*2 + gc.n_rays)
            self.rays = <gausslet_t*>realloc(self.rays, self.max_size*sizeof(gausslet_t))
        memcpy(self.rays + self.n_rays, gc.rays, gc.n_rays*sizeof(gausslet_t))
        self.n_rays += gc.n_rays
    
    @classmethod
    def from_array(cls, np_.ndarray data):
        """Creates a new RayCollection from the given numpy array. The array
        dtype should be a ctracer.ray_dtype. The data is copied into the 
        RayCollection
        """
        cdef: 
            int size=data.shape[0]
            GaussletCollection rc = GaussletCollection(size)
            
        if data.dtype != gausslet_dtype:
            raise ValueError("Array must have gausslet_dtype dtype")
        
        memcpy(rc.rays, <np_.float64_t *>data.data, size*sizeof(gausslet_t))
        rc.n_rays = size
        return rc
    
    @classmethod
    def from_rays(cls, np_.ndarray data):
        cdef: 
            int i, j, size=data.shape[0]
            GaussletCollection rc = GaussletCollection(size)
            gausslet_t *gc
            ray_t *ray
            para_t *p
            ray_t *_data = <ray_t *>data.data
            
        if data.dtype != ray_dtype:
            raise ValueError("Array must have gausslet_dtype dtype")
            
        for i in range(size):
            gc = rc.rays+i
            ray = &(_data[i])
            gc.base_ray = ray[0]
            for j in range(6):
                p = gc.para + j
                p.origin = ray.origin
                p.direction = ray.direction
                p.normal = ray.normal
                p.length = ray.length
        rc.n_rays = size
        return rc
    
    property lagrange_invariant:
        def __get__(self):
            cdef:
                int i,j
                unsigned long k, N = self.n_rays
                np_.ndarray out = np.empty((N,), dtype='d' )
                double[:] _out = out
                double v=0
                gausslet_t g
                ray_t r
                para_t p
                vector_t[6] h
                vector_t[6] u
                vector_t axis1, axis2, o
                double[:] wavelen = self.wavelengths
                
                
            for k in range(N):
                g = self.rays[k]
                r = g.base_ray
                axis1 = norm_(r.E_vector)
                axis2 = cross_(axis1, r.direction)
                
                for i in range(6):
                    p = g.para[i]
                    o = subvv_(p.origin, r.origin)
                    h[i].x = dotprod_(o, axis1)
                    h[i].y = dotprod_(o, axis2)
                    h[i].z = 0.0
                    u[i].x = dotprod_(p.direction, axis1)
                    u[i].y = dotprod_(p.direction, axis2)
                    u[i].z = 0.0
                    
                v = 0
                v += (dotprod_(h[0],u[5]) - dotprod_(h[5],u[0]))**2
                v += (dotprod_(h[1],u[2]) - dotprod_(h[2],u[1]))**2
                v += (dotprod_(h[3],u[4]) - dotprod_(h[4],u[3]))**2
                v += (dotprod_(h[1],u[4]) - dotprod_(h[4],u[1]))**2
                v += (dotprod_(h[0],u[3]) - dotprod_(h[3],u[0]))**2
                v += (dotprod_(h[2],u[5]) - dotprod_(h[5],u[2]))**2
                
                _out[k] = 1000*sqrt(v/6)/ wavelen[r.wavelength_idx]
            return out
        
    def project_to_plane(self, origin, direction):
        """
        Project the rays in the collection onto the intersection with the given plane,
        defined by an origin point on the plane and the plane normal vector.
        """
        cdef:
            vector_t o,d
            unsigned int i, j
            gausslet_t *gc
            ray_t *ray
            para_t *para
            double a
            #double complex[:] wl = self.wavelengths
            
        o.x, o.y, o.z = origin
        d.x, d.y, d.z = direction
        d = norm_(d)
        
        for i in range(self.n_rays):
            gc = &(self.rays[i])
            ray = &(gc.base_ray)
            a = dotprod_(subvv_(o, ray.origin),d) / dotprod_(ray.direction, d)
            ray.origin = addvv_(ray.origin, multvs_(ray.direction,a))
            ray.accumulated_path += ray.refractive_index.real * a
            ### Handle absorption along ray
            ### Am ignoring this for now. FIXME
            for j in range(6):
                para = &(gc.para[j])
                a = dotprod_(subvv_(o, para.origin),d) / dotprod_(para.direction, d)
                para.origin = addvv_(para.origin, multvs_(para.direction,a))
                
    def scale_amplitude(self, double complex scale):
        cdef:
            unsigned long i=self.n_rays
            
        for i in range(self.n_rays):
            self.rays[i].base_ray.E1_amp *= scale
            self.rays[i].base_ray.E2_amp *= scale
    
    def config_parabasal_rays(self, double[:] wavelength_list, double radius, double working_dist):
        """
        Initialise the parabasal rays for a symmetric (i.e. circular) modes, 
        using the base_ray data for wavelength, and the given beam waist 1/e^2 radius.
        'working_dist' indicates the distance from the base_ray origin to the centre 
        of the gaussian beam waist. Negative values imply a beam waist before the origin. 
        'radius' is given in mm.
        """
        cdef:
            int i,j
            gausslet_t *gc
            double theta0, denom, angle
            vector_t o, d, d1, d2, base_d, da, db
            
        for i in range(self.n_rays):
            gc = self.rays+i
            base_d = norm_(gc.base_ray.direction)
            if base_d.x > base_d.y:
                o.x=0.0
                o.y=1.0
                o.z=0.0
            else:
                o.x=1.0
                o.y=0.0
                o.z=0.0
            d1 = norm_(cross_(base_d, o))
            d2 = norm_(cross_(base_d, d1))
            ### Divergence of the gausslet
            theta0 = wavelength_list[gc.base_ray.wavelength_idx]/(M_PI*radius*1000.0)
            
            for j in range(0,6,2):
                angle = (j*2*M_PI/6)# + i*(2*M_PI)/self.n_rays
                o = addvv_(multvs_(d1, radius*cos(angle)),
                                multvs_(d2, radius*sin(angle)))
                o = addvv_(o, gc.base_ray.origin)
                o = addvv_(o, multvs_(base_d, working_dist))
                ###direction of ray
                d = addvv_(multvs_(d1, -theta0*sin(angle)),
                            multvs_(d2, theta0*cos(angle)))
                da = addvv_(base_d, d)
                db = subvv_(base_d, d)
                
                gc.para[j].direction = norm_(da)
                gc.para[j].origin = subvv_(o, multvs_(da,working_dist))
                
                gc.para[j+1].direction = norm_(db)
                gc.para[j+1].origin = subvv_(o, multvs_(db,working_dist))
                
                gc.para[j].normal = gc.base_ray.normal
                gc.para[j+1].normal = gc.base_ray.normal 
                gc.para[j].length = gc.base_ray.length
                gc.para[j+1].length = gc.base_ray.length
                
    property base_rays:
        def __get__(self):
            return GaussletBaseRayView(self)
        
    property total_power:
        def __get__(self):
            cdef:
                unsigned long i
                double pwr=0.0
                ray_t *ray
                double n
                
            for i in range(self.n_rays):
                ray = &(self.rays[i].base_ray)
                n = ray.refractive_index.real
                pwr += (ray.E1_amp.real * ray.E1_amp.real)*n
                pwr += (ray.E1_amp.imag * ray.E1_amp.imag)*n
                pwr += (ray.E2_amp.real * ray.E2_amp.real)*n
                pwr += (ray.E2_amp.imag * ray.E2_amp.imag)*n
            return pwr
                
    property wavelengths:
        def __get__(self):
            return np.asarray(self._wavelengths)
        
        def __set__(self, wl_list):
            self._wavelengths = np.asarray(wl_list, dtype=np.double)
        
    property para_origin:
        def __get__(self):
            cdef:
                np_.ndarray out = np.empty((self.n_rays,6,3), dtype='d')
                size_t i, j
                vector_t v
            for i in xrange(self.n_rays):
                for j in xrange(6):
                    v = self.rays[i].para[j].origin
                    out[i,j,0] = v.x
                    out[i,j,1] = v.y
                    out[i,j,2] = v.z
            return out
        
    property para_direction:
        def __get__(self):
            cdef:
                np_.ndarray out = np.empty((self.n_rays,6,3), dtype='d')
                size_t i, j
                vector_t v
            for i in xrange(self.n_rays):
                for j in range(6):
                    v = self.rays[i].para[j].direction
                    out[i,j,0] = v.x
                    out[i,j,1] = v.y
                    out[i,j,2] = v.z
            return out
        
    property para_normal:
        def __get__(self):
            cdef:
                np_.ndarray out = np.empty((self.n_rays,6,3), dtype='d')
                size_t i, j
                vector_t v
            for i in xrange(self.n_rays):
                for j in range(6):
                    v = self.rays[i].para[j].normal
                    out[i,j,0] = v.x
                    out[i,j,1] = v.y
                    out[i,j,2] = v.z
            return out
        
    property para_termination:
        def __get__(self):
            cdef:
                np_.ndarray aout = np.empty((self.n_rays,6,3), dtype='d')
                double [:,:,:] out = aout
                size_t i,j
                vector_t v
                para_t *p
            for i in xrange(self.n_rays):
                for j in range(1,7):
                    p = self.rays[i].para + j
                    v = addvv_(p.origin, 
                           multvs_(p.direction, 
                                    p.length))
                    out[i,j,0] = v.x
                    out[i,j,1] = v.y
                    out[i,j,2] = v.z
            return aout
        
    
cdef class InterfaceMaterial(object):
    """Abstract base class for objects describing
    the materials characterics of a Face
    """
    
    def __cinit__(self):
        self.wavelengths = np.array([], dtype=np.double)
    
    cdef void eval_child_ray_c(self, ray_t *old_ray, 
                                unsigned int ray_idx, 
                                vector_t p, 
                                orientation_t orient,
                                RayCollection new_rays):
        pass
    
    cdef para_t eval_parabasal_ray_c(self, ray_t *base_ray, 
                                     vector_t direction, #incoming ray direction
                                   vector_t point, #position of intercept
                                   orientation_t orient,
                                   unsigned int ray_type_id, #bool, if True, it's a reflected ray
                                   ):
        cdef:
            vector_t cosThetaNormal, reflected, normal
            para_t para_out
            double cosTheta
        
        normal = norm_(orient.normal)
        if ray_type_id & REFL_RAY:
            cosTheta = dotprod_(normal, direction)
            cosThetaNormal = multvs_(normal, cosTheta)
            reflected = subvv_(direction, multvs_(cosThetaNormal, 2))
            para_out.direction = reflected
        else:
            para_out.direction = direction
        para_out.origin = point
        para_out.normal = normal
        para_out.length = INF
        return para_out
    
    def is_decomp_material(self):
        return False
    
    cdef void eval_decomposed_rays_c(self, GaussletCollection child_rays):
        ### This function is called at the end of tracing to comput additional rays
        ### due to mode decomposition.
        return
    
    def eval_child_ray(self, Ray old_ray, ray_idx, point, 
                        normal, tangent, RayCollection new_rays):
        cdef:
            vector_t p
            orientation_t n
            Ray out=Ray()
            unsigned int idx
        
        p = set_v(point)
        n.normal = set_v(normal)
        n.tangent = set_v(tangent)
        self.eval_child_ray_c(&old_ray.ray, ray_idx, 
                                        p, n, new_rays)
        
    def eval_parabasal_ray(self, Ray base_ray, direction, point, normal, tangent, reflect=False):
        cdef:
            vector_t d=set_v(direction), p = set_v(point)
            orientation_t n
            ParabasalRay out = ParabasalRay()
            unsigned int ray_type=1 if reflect else 0
            
        n.normal = set_v(normal)
        n.tangent = set_v(tangent)
        out.ray = self.eval_parabasal_ray_c(&base_ray.ray, d, p, n, ray_type)
        return out
            
        
    property wavelengths:
        def __set__(self, double[:] wavelengths):
            self._wavelengths = wavelengths
            self.on_set_wavelengths()
            
        def __get__(self):
            return self._wavelengths
        
    cdef on_set_wavelengths(self):
        pass
    
    
cdef class Shape:        
    cdef bint point_inside_c(self, double x, double y):
        return 1
    
    def point_inside(self, double x, double y):
        return self.point_inside_c(x,y) 
    
    
cdef class Distortion:
    """A abstract base class to represents distortions on a face, a z-offset 
    as a function of (x,y).
    """
    cdef vector_t z_offset_and_gradient_c(self, double x, double y) nogil:
        """The z-axis surface sag is returned as the z-component 
        of the output vector. The x- and y-components of the surface
        gradient are placed in the x- and y- components of vector.
        I.e. the returned vector is
            (dz(x,y)/dx, dz(x,y)/dy, z(x,y))
        """
        cdef:
            vector_t p
        p.x=0.0
        p.y=0.0
        p.z=0.0
        return p
    
    cdef double z_offset_c(self, double x, double y) nogil:
        return 0.0
    
    def z_offset_and_gradient(self, double[:] x, double[:] y):
        cdef:
            vector_t v
            size_t i, n=len(x)
            double[:,:] out = np.empty((n, 3), 'd')
            
        if len(y) != n:
            raise ValueError("Both x and y must have the same length")
            
        for i in range(n):
            v = self.z_offset_and_gradient_c(x[i],y[i])
            out[i,0] = v.x
            out[i,1] = v.y
            out[i,2] = v.z
        return np.asarray(out)
    
    def z_offset(self, double[:] x, double[:] y):
        cdef:
            size_t i, n=len(x)
            double[:] out = np.empty((n,), 'd')
            
        if len(y) != n:
            raise ValueError("Both x and y must have the same length")
            
        for i in range(n):
            out[i] = self.z_offset_c(x[i],y[i])
        return np.asarray(out)
    
    
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
            from .cmaterials import PECMaterial
            self.material = PECMaterial()
        self.invert_normal = int(kwds.get('invert_normal', 0))
        
    
    cdef double intersect_c(self, vector_t p1, vector_t p2, int is_base_ray):
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
    
    def intersect(self, p1, p2, int is_base_ray):
        cdef:
            vector_t p1_, p2_
            double dist
        
        p1_ = set_v(p1)
        p2_ = set_v(p2)
        dist = self.intersect_c(p1_, p2_, is_base_ray)
        return dist

    cdef vector_t compute_normal_c(self, vector_t p):
        return p
    
    cdef vector_t compute_tangent_c(self, vector_t p):
        cdef vector_t tangent
        tangent.x = 1.0
        tangent.y = 0.0
        tangent.z = 0.0
        return tangent
    
    def compute_normal(self, p):
        """Compute normal vector at a given point, in local
        face coordinates
        """
        cdef vector_t p_, n
        p_.x, p_.y, p_.z = p
        n = self.compute_normal_c(p_)
        return (n.x, n.y, n.z)
    
    def compute_tangent(self, p):
        """Compute the surface tangent at a given point,
        in face-local coordinates"""
        cdef vector_t tanget, p_
        
        p_.x, p_.y, p_.z = p
        tangent = self.compute_tangent_c(p_)
        return (tangent.x, tangent.y, tangent.z)
        


cdef class FaceList(object):
    """A group of faces which share a transform"""
    def __cinit__(self, owner=None):
        self.transform = Transform()
        self.inverse_transform = Transform()
        self.owner = owner
        
    cpdef void sync_transforms(self):
        """sets the transforms from the owner's VTKTransform
        """
        try:
            trans = self.owner.transform
        except AttributeError:
            print("NO OWNER", self.owner)
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
        
     
    cdef int intersect_c(self, ray_t *ray, vector_t ray_end):
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
            dist = face.intersect_c(p1, p2, 1)
            if face.tolerance < dist < ray.length:
                ray.length = dist
                all_idx = face.idx
                ray.end_face_idx = all_idx
        return all_idx
    
    def intersect(self, Ray r):
        cdef vector_t P1_
        cdef int idx
        
        P1_ = addvv_(r.ray.origin, multvs_(r.ray.direction, r.ray.length))
        idx = self.intersect_c(&r.ray, P1_)
        return idx
    
    cdef int intersect_para_c(self, para_t *ray, vector_t ray_end, Face face):
        cdef:
            vector_t p1 = transform_c(self.inv_trans, ray.origin)
            vector_t p2 = transform_c(self.inv_trans, ray_end)
            unsigned int i
            double dist
        
        dist = face.intersect_c(p1, p2, 0)
        #print("p1:", p1.x, p1.y, p1.z, "p2:", p2.x, p2.y, p2.z, "rlen:", ray.length, "dist:", dist)
        if face.tolerance < dist < ray.length:
            ray.length = dist
            return 0
        return -1
    
    def intersect_para(self, ParabasalRay r, Face face):
        cdef vector_t P1_
        cdef int idx
        
        P1_ = addvv_(r.ray.origin, multvs_(r.ray.direction, r.ray.length))
        idx = self.intersect_para_c(&r.ray, P1_, face)
        return idx
    
    cdef orientation_t compute_orientation_c(self, Face face, vector_t point):
        cdef orientation_t out
        
        point = transform_c(self.inv_trans, point)
        out.normal = face.compute_normal_c(point)
        out.tangent = face.compute_tangent_c(point)
        if face.invert_normal:
            out.normal = invert_(out.normal)
            out.tangent = invert_(out.tangent)
        out.normal = rotate_c(self.trans, out.normal)
        out.tangent = rotate_c(self.trans, out.tangent)
        return out
    
    def compute_orientation(self, Face face, point):
        cdef:
            vector_t p
            orientation_t o
        p = set_v(point)
        o = self.compute_orientation_c(face, p)
        return (o.normal.x, o.normal.y, o.normal.z), (o.tangent.x, o.tangent.y, o.tangent.z)
    

##################################
### Python module functions
##################################

cdef double ray_power_(ray_t ray):
    cdef double P1, P2
    
    P1 = (ray.E1_amp.real**2 + ray.E1_amp.imag**2)*ray.refractive_index.real
    P2 = (ray.E2_amp.real**2 + ray.E2_amp.imag**2)*ray.refractive_index.real
    ### We don't need to handle incident aspect area here as the ray already compensate for this
    #aspect = abs(dotprod(norm(ray.direction), normal))
    return (P1+P2) 
            
            
def select_ray_intersections(FaceList face_set, list ray_col_list):
    cdef:
        size_t i,j
        ray_t ray
        int idx, wl_offset=0
        vector_t point
        RayCollection rc, rc_out
        list wl_list=[]
        np_.int64_t[:] inverse
        
    rc_out = RayCollection(len(ray_col_list[0]))
        
    for j in range(len(ray_col_list)):
        rc = <RayCollection>ray_col_list[j]
        wl_list.append(rc.wavelengths)
        
        for i in range(rc.n_rays):
            ray = rc.rays[i] #need a copy of the ray to prevent mutating it in place
            point = addvv_(ray.origin, 
                                multvs_(ray.direction, 
                                        ray.length))
            
            idx = (<FaceList>face_set).intersect_c(&ray, point)
            if idx >= 0:
                ray.wavelength_idx += wl_offset
                rc_out.add_ray_c(ray)
        wl_offset += len(rc.wavelengths)
        
    reduced, inverse = np.unique(np.concatenate(wl_list), return_inverse=True)
    for i in range(len(rc_out)):
        idx = rc_out.rays[i].wavelength_idx
        rc_out.rays[i].wavelength_idx = inverse[idx]
        
    rc_out.wavelengths = reduced
    return rc_out
    
    
def select_gausslet_intersections(FaceList face_set, list ray_col_list):
    cdef:
        size_t i,j
        gausslet_t g
        ray_t *ray
        int idx, wl_offset=0
        vector_t point
        GaussletCollection rc, rc_out
        list wl_list=[]
        np_.int64_t[:] inverse
        
    rc_out = GaussletCollection(len(ray_col_list[0]))
        
    for j in range(len(ray_col_list)):
        rc = <GaussletCollection>ray_col_list[j]
        wl_list.append(rc.wavelengths)
        
        for i in range(rc.n_rays):
            g = rc.rays[i] #need a copy of the ray to prevent mutating it in place
            ray = &g.base_ray
            point = addvv_(ray.origin, 
                                multvs_(ray.direction, 
                                        ray.length))
            
            idx = (<FaceList>face_set).intersect_c(ray, point)
            if idx >= 0:
                ray.wavelength_idx += wl_offset
                rc_out.add_gausslet_c(g)
        wl_offset += len(rc.wavelengths)
        
    reduced, inverse = np.unique(np.concatenate(wl_list), return_inverse=True)
    for i in range(len(rc_out)):
        idx = rc_out.rays[i].base_ray.wavelength_idx
        rc_out.rays[i].base_ray.wavelength_idx = inverse[idx]
        
    rc_out.wavelengths = reduced
    return rc_out
    


cdef RayCollection trace_segment_c(RayCollection rays, 
                                    list face_sets, 
                                    list all_faces,
                                    list decomp_faces,
                                    float max_length):
    cdef:
        FaceList face_set #a FaceList
        Face face
        size_t i, j, n_sets=len(face_sets)
        vector_t point
        orientation_t orient
        int idx, nearest_set=-1, nearest_idx=-1
        ray_t *ray
        RayCollection new_rays
   
    #need to allocate the output rays here 
    new_rays = RayCollection(rays.n_rays)
    new_rays.parent = rays
    
    for i in range(rays.n_rays):
        ray = rays.rays + i
        ray.length = max_length
        ray.end_face_idx = -1
        nearest_idx=-1
        point = addvv_(ray.origin, 
                            multvs_(ray.direction, 
                                    max_length))
        #print "points", P1, P2
        for j in range(n_sets):
            face_set = face_sets[j]
            #intersect_c returns the face idx of the intersection, or -1 otherwise
            idx = (<FaceList>face_set).intersect_c(ray, point)
            if idx >= 0:
                nearest_set = j
                nearest_idx = idx
        if nearest_idx >= 0:
            #print "GET FACE", nearest.face_idx, len(all_faces)
            face = all_faces[nearest_idx]
            face.count += 1
            #print "ray length", ray.length
            point = addvv_(ray.origin, multvs_(ray.direction, ray.length))
            orient = (<FaceList>(face_sets[nearest_set])).compute_orientation_c(face, point)
            #print "s normal", normal
            (<InterfaceMaterial>(face.material)).eval_child_ray_c(ray, i, 
                                                    point,
                                                    orient,
                                                    new_rays
                                                    )
    return new_rays


def trace_segment(RayCollection rays, 
                    list face_sets, 
                    list all_faces,
                    max_length=100,
                    decomp_faces=[]):
    for fs in face_sets:
        fs.sync_transforms()
    return trace_segment_c(rays, face_sets, all_faces, decomp_faces, max_length)

def trace_gausslet(GaussletCollection rays, 
                    list face_sets, 
                    list all_faces,
                    max_length=100,
                    decomp_faces=[]):
    for fs in face_sets:
        (<FaceList>fs).sync_transforms()
    return trace_gausslet_c(rays, face_sets, all_faces, decomp_faces, max_length)


cdef GaussletCollection trace_gausslet_c(GaussletCollection gausslets, 
                                    list face_sets, 
                                    list all_faces,
                                    list decomp_faces,
                                    double max_length):
    cdef:
        FaceList face_set #a FaceList
        Face face
        size_t i, j, n_sets=len(face_sets), n_decomp = len(decomp_faces)
        vector_t point
        orientation_t orient
        int idx, nearest_set=-1, nearest_idx=-1
        gausslet_t *gausslet
        ray_t *ray
        GaussletCollection new_gausslets
        RayCollection child_rays
   
    #need to allocate the output rays here 
    new_gausslets = GaussletCollection(gausslets.n_rays)
    new_gausslets.parent = gausslets
    
    child_rays = RayCollection(2)
    
    for i in range(gausslets.n_rays):
        gausslet = gausslets.rays + i
        ray = &gausslet.base_ray
        ray.end_face_idx = -1
        nearest_idx=-1
        point = addvv_(ray.origin, 
                            multvs_(ray.direction, 
                                    max_length))
        #print "points", P1, P2
        for j in range(n_sets):
            face_set = face_sets[j]
            #intersect_c returns the face idx of the intersection, or -1 otherwise
            idx = (<FaceList>face_set).intersect_c(ray, point)
            if idx >= 0:
                nearest_set = j
                nearest_idx = idx
        if nearest_idx >= 0:
            #print "GET FACE", nearest.face_idx, len(all_faces)
            face = all_faces[nearest_idx]
            face.count += 1
            #print "ray length", ray.length
            point = addvv_(ray.origin, multvs_(ray.direction, ray.length))
            face_set = (<FaceList>(face_sets[nearest_set])) 
            orient = face_set.compute_orientation_c(face, point)
            #print "s normal", normal
            ### Clear the child_rays structure
            child_rays.n_rays = 0
            (<InterfaceMaterial>(face.material)).eval_child_ray_c(ray, i, 
                                                    point,
                                                    orient,
                                                    child_rays
                                                    )
            trace_parabasal_rays(gausslet, child_rays, face, face_set, new_gausslets, max_length)
            
    for j in range(n_decomp):
        face = decomp_faces[j]
        if face.count:
            (<InterfaceMaterial>(face.material)).eval_decomposed_rays_c(new_gausslets)
            face.count = 0
            
    new_gausslets.reset_length_c(max_length)
    return new_gausslets


cdef void trace_parabasal_rays(gausslet_t *g_in, RayCollection base_rays, Face face, FaceList face_set, 
                          GaussletCollection new_gausslets, double max_length):
    cdef:
        unsigned int i,j
        para_t *para_ray
        gausslet_t gausslet
        vector_t ray_end, point[6]
        orientation_t orient[6]
        InterfaceMaterial material = (<InterfaceMaterial>(face.material))
        
    for j in range(6):
        para_ray = g_in.para + j
        ray_end = addvv_(para_ray.origin, 
                        multvs_(para_ray.direction, 
                                max_length))
        if face_set.intersect_para_c(para_ray, ray_end, face):
            return
        point[j] = addvv_(para_ray.origin, multvs_(para_ray.direction, para_ray.length))
        orient[j] = face_set.compute_orientation_c(face, point[j])
    
    for i in range(base_rays.n_rays):
        gausslet.base_ray = base_rays.rays[i]
        for j in range(6):
            para_ray = g_in.para + j            
            gausslet.para[j] = material.eval_parabasal_ray_c(base_rays.rays + i,
                                                            para_ray.direction, #incoming ray direction
                                                            point[j], #position of intercept
                                                            orient[j],
                                                            gausslet.base_ray.ray_type_id, #indicates if it's a transmitted or reflected ray 
                                                            )
        new_gausslets.add_gausslet_c(gausslet)
            
            

def transform(Transform t, p):
    cdef vector_t p1, p2
    assert isinstance(t, Transform)
    assert len(p)==3
    p1.x = p[0]
    p1.y = p[1]
    p1.z = p[2]
    p2 = transform_c(t.trans, p1)
    return (p2.x, p2.y, p2.z)
    
def get_ray_size():
    return sizeof( ray_t )

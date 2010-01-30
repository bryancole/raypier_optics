#!/usr/bin/env python

cdef extern from "math.h":
    double sqrt
    double INFINITY
    
ctypedef double[3] vector

cdef double sep(vector p1, vector p2):
    cdef double val
    val = 0
    for i in range(3):
        val += (p2[i] - p1[i])**2
    return sqrt(val)


cdef struct Intersection:
    Face face #the intersecting face
    vector point #position of intersection
    double dist #fractional of ray between origin and intersection
    
cdef struct complex:
    double real
    double imag


cdef class Ray:
    #vectors
    cdef vector origin, direction, normals, E_vector
    #complex attribs
    cdef complex refractive_index, E1_amp, E2_amp
    #simple attribs
    cdef double length, wavelength
    #objects
    cdef object face, end_face, child_refl, child_trans
    
    
cdef class Face(object):
    cdef public object owner
    cdef public char *name
    cdef public double tolerance
    
    def __cinit__(self):
        self.name = "base Face class"
        self.tolerance = 0.0001
    
    cdef vector intersect_c(self, vector p1, vector p2):
        """returns the intersection in terms of the 
        fractional distance between p1 and p2.
        p1 and p2 are in the local coordinate system
        """
        return 0.5

    cdef vector compute_normal_c(self, vector p):
        return Vector(p.z,p.y,p.x)
    
    cdef Ray eval_child_ray_c(self, Ray ray, vector p):
        return ray


cdef class FaceList(object):
    """A group of faces which share a transform"""
    cdef public object transform
    cdef public object inverse_transform
    cdef public list faces
     
    cdef Intersection intersect_c(self, vector P1, vector P2, double max_length):
        """Finds the face with the nearest intersection
        point, for the ray defined by the two input points,
        P1 and P2 (in global coords).
        """
        cdef vector p1, p2, point
        cdef list faces
        cdef double d, dist=INFINITY
        cdef Face nearest=None
        cdef Intersection inter
        
        p1 = transform(self.transform, P1)
        p2 = transform(self.transform, P2)
        
        faces = self.faces
        
        for i in xrange(len(faces)):
            point = (<Face>(faces[i])).intersect_c(p1, p2)
            d = sep(p1, point)
            if 0.0 < d < dist:
                dist = d
                nearest = <Face>(faces[i])
        
        inter.face = nearest
        inter.point = transform(self.inverse_transform, point)
        inter.dist = dist/sep(p1,p2)
        return inter


cdef Ray trace_ray(Ray ray, list face_sets, double max_length):
    cdef vector P1, P2, direction, normal
    cdef double dist=INFINITY
    cdef Intersection inter, nearest
    cdef unsigned int i, n_face_sets=len(face_sets)

    P1 = ray.origin
    direction = ray.direction
    for i in xrange(3):
        P2[i] = P1[i] + (direction[i]*max_length)

    for i in xrange(n_face_sets): #need to make this a C array
        inter = (<FaceList>face_sets[i]).intersect_c(P1, P2, max_length)
        if 0.0 < inter.dist < dist:
            nearest = inter
    
    #evaluate new ray
    return (<Face>(inter.face)).eval_new_ray_c(ray, nearest.point)






cdef object trace_segment(object rays, list optics):
    FaceList face_set
    int size, i
    float[3] P1, P2
    float max_length
    Intersection inter, inter2
    
    face_sets = [o.faces for o in optics]
    size = rays.number
    origin = rays.origin
    direction = rays.direction
    max_length = rays.max_length
    
    #need to allocate the output rays here
    
    for i in range(size):
        P1 = origin[i]
        P2 = P1 + (direction[i] * max_length)
        inter = (<FaceList>(face_sets[0])).intersect_c(P1, P2)
        dist = sep(inter.point, P1)
        for face_set in face_sets[1:]: #need to make this a C array
            inter2 = (<FaceList>face_set).intersect_c(P1, P2)
            dist2 = sep(inter2.point, P1)
            if dist2 < dist:
                dist = dist2
                inter = inter2
        
        #now compute normal
        normal = inter.face.compute_normal(inter.point)
        
        #evaluate new ray
        
                
        
        
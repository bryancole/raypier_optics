
cdef extern from "math.h":
    double sqrt(double arg)
    double fabs(double arg)
    double INFINITY

from ctracer cimport InterfaceMaterial, norm_, dotprod_, \
        multvs_, subvv_, vector_t, ray_t, RayCollection, \
        complex_t, mag_sq_, Ray, cross_, set_v


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
    n = set_v(normal)
    r = convert_to_sp(ray.ray, n)
    out = Ray()
    out.ray = r
    return out


cdef class OpaqueMaterial(InterfaceMaterial):
    """A perfect absorber i.e. it generates no rays
    """
    cdef eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            vector_t normal,
                            RayCollection new_rays):
        pass
    
    
cdef class TransparentMaterial(InterfaceMaterial):
    """A perfect transmisser i.e. it generates an
    outgoing ray with identical direction, polarisation etc.
    to the incoming ray. It does project the polarisation
    vectors to it's S- and P-directions, however.
    """
    cdef eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            vector_t normal,
                            RayCollection new_rays):
        cdef:
            vector_t cosThetaNormal, reflected
            ray_t sp_ray
            complex_t cpx
            double cosTheta
        
        normal = norm_(normal)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.parent_idx = idx
        new_rays.add_ray_c(sp_ray)


cdef class PECMaterial(InterfaceMaterial):
    """Simulates a Perfect Electrical Conductor
    """
    cdef eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            vector_t normal,
                            RayCollection new_rays):
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
        sp_ray = convert_to_sp(in_ray[0], normal)
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
        new_rays.add_ray_c(sp_ray)
    
    
cdef class DielectricMaterial(InterfaceMaterial):
    """Simulates Fresnel reflection and refraction at a
    normal dielectric interface
    """
    cdef complex_t n_inside_, n_outside_
    
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
            
    cdef eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            vector_t normal,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2, in_direction
            ray_t sp_ray
            complex_t cpx
            double cosTheta, n1, n2, N2, cos1
            double N2cosTheta, N2_sin2, tan_mag_sq, c2
            double cos2, Two_n1_cos1, aspect, T_p, T_s
            int flip
            
        normal = norm_(normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray[0], normal)
        cosTheta = dotprod_(normal, in_direction)
        cos1 = fabs(cosTheta)
        
        #print "TRACE"
        #print normal, in_direction
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside_.real
            n2 = self.n_inside_.real
            sp_ray.refractive_index = self.n_inside_
            flip = 1
            #print "out to in", n1, n2
        else:
            n1 = self.n_inside_.real
            n2 = self.n_outside_.real
            sp_ray.refractive_index = self.n_outside_
            flip = -1
            #print "in to out", n1, n2
            
        N2 = (n2/n1)**2
        N2cosTheta = N2*cos1
        
        N2_sin2 = (cosTheta*cosTheta) + (N2 - 1)
        #print "TIR", N2_sin2, cosTheta, N2, cos1
        #print (normal.x, normal.y, normal.z), in_direction
        cosThetaNormal = multvs_(normal, cosTheta)
        if N2_sin2 < 0.0:
            #total internal reflection
            reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INFINITY
            sp_ray.E1_amp.real *= -1
            sp_ray.E1_amp.imag *= -1
            sp_ray.E2_amp.real *= -1
            sp_ray.E2_amp.imag *= -1
            sp_ray.parent_idx = idx
        else:
            #normal transmission            
            tangent = subvv_(in_direction, cosThetaNormal)
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
        new_rays.add_ray_c(sp_ray)
        
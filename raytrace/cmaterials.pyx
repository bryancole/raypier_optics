
cdef extern from "math.h":
    double sqrt(double arg)
    double fabs(double arg)
    double sin(double arg)
    double cos(double arg)
    double atan2(double y, double x)
    double M_PI
    
cdef extern from "float.h":
    double DBL_MAX
    
cdef extern from "complex.h":
    double complex csqrt(double complex)
    double cabs(double complex)
    
cdef double INF=(DBL_MAX+DBL_MAX)

from ctracer cimport InterfaceMaterial, norm_, dotprod_, \
        multvs_, subvv_, vector_t, ray_t, RayCollection, \
        complex_t, mag_sq_, Ray, cross_, set_v, ray_power_


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
    E2_vector = norm_(cross_(ray.direction, ray.E_vector))
    E1_vector = norm_(cross_(E2_vector, ray.direction))
    normal = norm_(normal)
        
    S_vector = cross_(ray.direction, normal)
    if S_vector.x==0 and S_vector.y==0 and S_vector.z==0:
        print "Normal incidence!"
        return ray
    print "Not normal incidence", S_vector.x, S_vector.y, S_vector.z
    S_vector = norm_(S_vector)
        
    v = cross_(ray.direction, S_vector)
    P_vector = norm_(v)
    
    A = dotprod_(E1_vector, S_vector)
    B = dotprod_(E2_vector, S_vector)
    
    S_amp.real = E1_amp.real*A + E2_amp.real*B
    S_amp.imag = E1_amp.imag*A + E2_amp.imag*B
    
    B = dotprod_(E1_vector, P_vector)
    A = dotprod_(E2_vector, P_vector)

    P_amp.real = E1_amp.real*B + E2_amp.real*A
    P_amp.imag = E1_amp.imag*B + E2_amp.imag*A
    
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
    """A perfect transmitter i.e. it generates an
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
        #sp_ray.E2_amp.real = -sp_ray.E2_amp.real
        #sp_ray.E2_amp.imag = -sp_ray.E2_amp.imag
        sp_ray.parent_idx = idx
        new_rays.add_ray_c(sp_ray)
        
        
cdef class LinearPolarisingMaterial(InterfaceMaterial):
    """Simulates a perfect polarising beam splitter. P-polarisation
    is 100% transmitted while S- is reflected"""
    
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
            vector_t tangent, tg2, in_direction
            ray_t sp_ray, sp_ray2
            double cosTheta
            
        normal = norm_(normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = sp_ray2 = convert_to_sp(in_ray[0], normal)
        cosTheta = dotprod_(normal, in_direction)            
        cosThetaNormal = multvs_(normal, cosTheta)
        
        #Reflect the S-polarisation
        reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.direction = reflected
        sp_ray.length = INF
        #set P-polarisation to zero
        sp_ray.E2_amp.real = 0.0
        sp_ray.E2_amp.imag = 0.0
        sp_ray.parent_idx = idx
        new_rays.add_ray_c(sp_ray)
            
        #Transmit the P-polarisation
        sp_ray2.origin = point
        sp_ray2.normal = normal
        sp_ray2.direction = in_direction
        sp_ray2.length = INF
        #set S-polarisation to zero
        sp_ray2.E1_amp.real = 0.0
        sp_ray2.E1_amp.imag = 0.0
        sp_ray2.parent_idx = idx
        new_rays.add_ray_c(sp_ray2)
        
        
cdef class WaveplateMaterial(InterfaceMaterial):
    """An idealised optical retarder"""
    cdef complex_t retardance_
    cdef vector_t fast_axis_
    
    property retardance:
        def __get__(self):
            cdef double val
            val = atan2(self.retardance_.imag, 
                         self.retardance_.real)/(2*M_PI)
            if val < 0:
                val += 1.0
            return val
            
        def __set__(self, val):
            self.retardance_.real = cos(val*2*M_PI)
            self.retardance_.imag = sin(val*2*M_PI)
    
    property fast_axis:
        def __get__(self):
            cdef vector_t ax = self.fast_axis_
            return (ax.x, ax.y, ax.z)
        
        def __set__(self, ax):
            self.fast_axis_.x = ax[0]
            self.fast_axis_.y = ax[1]
            self.fast_axis_.z = ax[2]
            
            
    def __cinit__(self, **kwds):
        self.retardance = kwds.get("retardance", 0.25)
        self.fast_axis = kwds.get("fast_axis", (1.0,0,0))
            
            
    cdef ray_t apply_retardance_c(self, ray_t r):
        cdef:
            complex_t E1=r.E1_amp, retard=self.retardance_
        r.E1_amp.real = E1.real*retard.real - E1.imag*retard.imag
        r.E1_amp.imag = E1.imag*retard.real + E1.real*retard.imag
        #Leave E2_amp untouched.
        return r
    
    def apply_retardance(self, Ray r):
        cdef Ray out
        out = Ray()
        out.ray = self.apply_retardance_c(r.ray)
        return out
            
    
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
            ray_t out_ray
            complex_t E1, retard
            
        normal = norm_(normal)
        in_direction = norm_(in_ray.direction)
        ###The "P" output of convert_to_sp with be aligned with the fast_axis
        ###The "S" output will thus be orthogonal to the fast axis
        out_ray = convert_to_sp(in_ray[0], self.fast_axis_)
        
        out_ray.origin = point
        out_ray.normal = normal
        out_ray.direction = in_direction
        out_ray.length = INF
        out_ray.parent_idx = idx
        
        out_ray = self.apply_retardance_c(out_ray)
        
        new_rays.add_ray_c(out_ray)
    
    
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
            sp_ray.length = INF
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
            sp_ray.length = INF
            sp_ray.E1_amp.real *= T_s
            sp_ray.E1_amp.imag *= T_s
            sp_ray.E2_amp.real *= T_p
            sp_ray.E2_amp.imag *= T_p
            sp_ray.parent_idx = idx
        new_rays.add_ray_c(sp_ray)
        

cdef class FullDielectricMaterial(DielectricMaterial):
    """Model for dielectric using full Fresnel equations 
    to give true phase and amplitude response
    """
    cdef public double reflection_threshold, transmission_threshold
    
    def __cinit__(self, **kwds):
        self.reflection_threshold = kwds.get('reflection_threshold', 0.1)
        self.transmission_threshold = kwds.get('transmission_threshold', 0.1)
    
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
            double complex n1, n2, cos2, sin2, R_p, R_s, T_p, T_s, E1_amp, E2_amp
            double cosTheta, cos1, P_in
            double tan_mag_sq, c2
            double Two_n1_cos1, aspect
            int flip
            
        normal = norm_(normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray[0], normal)
        E1_amp = sp_ray.E1_amp.real + 1.0j*sp_ray.E1_amp.imag
        E2_amp = sp_ray.E2_amp.real + 1.0j*sp_ray.E2_amp.imag
        cosTheta = dotprod_(normal, in_direction)
        cos1 = fabs(cosTheta)
        sin1 = sqrt(1 - cos1*cos1)
        
        #print "TRACE"
        #print normal, in_direction
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside_.real + 1.0j*self.n_outside_.imag
            n2 = self.n_inside_.real + 1.0j*self.n_inside_.imag
            flip = 1
            #print "out to in", n1, n2
        else:
            n1 = self.n_inside_.real + 1.0j*self.n_inside_.imag
            n2 = self.n_outside_.real + 1.0j*self.n_outside_.imag
            flip = -1
            #print "in to out", n1, n2
            
        #apply Snell's law. These become complex.
        sin2 = (n1*sin1/n2)
        cos2 = csqrt(1 - sin2*sin2)
        #print "cos1", cos1
        #print "cos2", cos2
        
        #print "TIR", N2_sin2, cosTheta, N2, cos1
        #print (normal.x, normal.y, normal.z), in_direction
        cosThetaNormal = multvs_(normal, cosTheta)
        
        #incoming power
        P_in =  n1.real*(E1_amp.real**2 + E1_amp.imag**2 + \
                         E2_amp.real**2 + E2_amp.imag**2)
            
        #Fresnel equations for reflection
        R_p = -(n2*cos1 - n1*cos2)/(n2*cos1 + n1*cos2)
        R_s = -(n2*cos2 - n1*cos1)/(n2*cos2 + n1*cos1)
        
        #modify in place to get reflected amplitudes
        R_s *= E1_amp
        R_p *= E2_amp
        
        #print "P_in:", P_in
        #print "P_R:", (cabs(R_s)**2 + cabs(R_p)**2)
        
        if ( n1.real*(cabs(R_s)**2 + cabs(R_p)**2)/P_in ) > self.reflection_threshold:
            reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INF
            sp_ray.E1_amp.real = R_s.real
            sp_ray.E1_amp.imag = R_s.imag
            sp_ray.E2_amp.real = R_p.real
            sp_ray.E2_amp.imag = R_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n1.real
            sp_ray.refractive_index.imag = n1.imag
            new_rays.add_ray_c(sp_ray)
            
        #normal transmission            
        tangent = subvv_(in_direction, cosThetaNormal)
        tg2 = multvs_(tangent, n1.real/n2.real) #This is an approximation for complex N
        tan_mag_sq = mag_sq_(tg2)
        c2 = sqrt(1-tan_mag_sq)
        transmitted = subvv_(tg2, multvs_(normal, c2*flip))
        
        aspect = sqrt(cos2.real/cos1)
        
        #Fresnel equations for transmission
        T_p = aspect * (2.0*cos1*n1) / ( n2*cos1 + n1*cos2 )
        T_s = aspect * (2.0*cos1*n1) / ( n2*cos2 + n1*cos1 )
        #print "T_s", T_s, "T_p", T_p
        
        #modify in place to get reflected amplitudes
        T_s *= E1_amp
        T_p *= E2_amp
        
        if ( n2.real*(cabs(T_s)**2 + cabs(T_p)**2)/P_in ) > self.transmission_threshold:
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = transmitted
            sp_ray.length = INF
            sp_ray.E1_amp.real = T_s.real
            sp_ray.E1_amp.imag = T_s.imag
            sp_ray.E2_amp.real = T_p.real
            sp_ray.E2_amp.imag = T_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n2.real
            sp_ray.refractive_index.imag = n2.imag
            new_rays.add_ray_c(sp_ray)
# cython: cdivision=True
cimport cython

cdef extern from "math.h":
    double sqrt(double arg)
    double fabs(double arg)
    double sin(double arg)
    double cos(double arg)
    double atan2(double y, double x)
    double erf(double arg)
    double M_PI
    
cdef extern from "float.h":
    double DBL_MAX

IF UNAME_SYSNAME == "Windows":
    cdef extern from "complex.h":
        double complex csqrt "sqrt" (double complex)
        double cabs "abs" (double complex)
        double complex cexp "exp" (double complex)
	
    cdef double complex I = 1j
ELSE:
    cdef extern from "complex.h":
        double complex csqrt (double complex)
        double cabs (double complex)
        double complex cexp (double complex)
        double complex I
    
cdef:
    double INF=(DBL_MAX+DBL_MAX)
    double SP_TOL=1.0e-10

from .ctracer cimport InterfaceMaterial, norm_, dotprod_, \
        multvs_, subvv_, vector_t, ray_t, RayCollection, \
        complex_t, mag_sq_, Ray, cross_, set_v, ray_power_, \
        orientation_t, addvv_, REFL_RAY, PARABASAL, \
        para_t, gausslet_t, GaussletCollection

import numpy as np
cimport numpy as np_


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
    #print("E_in:", E1_amp.real, E1_amp.imag, E2_amp.real, E2_amp.imag)
    #print("E_vector:", E1_vector.x, E1_vector.y, E1_vector.z, E2_vector.x, E2_vector.y, E2_vector.z)
    S_vector = cross_(ray.direction, normal)
    if fabs(S_vector.x)<SP_TOL and fabs(S_vector.y)<SP_TOL and fabs(S_vector.z)<SP_TOL:
        return ray
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
    
    #ray.phase = 0.0
    return ray


def Convert_to_SP(Ray ray, normal):
    cdef vector_t n
    cdef ray_t r
    n = set_v(normal)
    r = convert_to_sp(ray.ray, n)
    out = Ray()
    out.ray = r
    return out


ctypedef double (*dispersion_curve) (double, double[:])


cdef double nondispersive_0(double wavelen, double[:] coefs):
    return coefs[0]

cdef double sellmeier_1(double wavelen, double[:] coefs ):
    cdef:
        double n2=1.0
        double wl2=wavelen*wavelen
        int i, n_coefs = coefs.shape[0]
        
    n2 += coefs[0]
    for i in range((n_coefs-1)/2):
        n2 += coefs[2*i + 1] * wl2 / (wl2 - coefs[2*i + 2]**2)
    return sqrt(n2)

cdef double sellmeier_2(double wavelen, double[:] coefs ):
    cdef:
        double n2=1.0
        double wl2=wavelen*wavelen
        int i, n_coefs = coefs.shape[0]
        
    n2 += coefs[0]
    for i in range((n_coefs-1)/2):
        n2 += coefs[2*i + 1] * wl2 / (wl2 - coefs[2*i + 2])
    return sqrt(n2)

### Polynomia formula
cdef double sellmeier_3(double wavelen, double[:] coefs ):
    cdef:
        double n2=coefs[0]
        int i, n_coefs = coefs.shape[0]
        
    for i in range(1,n_coefs-1,2):
        n2 += coefs[i]*( wavelen**(coefs[i+1]) )
    return sqrt(n2)
        

cdef double sellmeier_5(double wavelen, double[:] coefs):
    cdef:
        double n= coefs[0]
        int i
        
    for i in range((coefs.shape[0]-1)/2):
        n += coefs[2*i + 1] * (wavelen**coefs[2*i+2])
        
        
cdef class BaseDispersionCurve(object):
    """
    Base class for DispersionCurve objects. This extension class provides 
    efficient C-API methods for evaluation of the refractive index for a
    given wavelength, based on one of various dispersion functions (e.g. 
    Sellmeier curves) and a set of coefficients obtained from https://refractiveindex.info
    
    :param formula_id: Specifies the formula to use in evaluating this dispersion-curve
    :type formula_id: int, readonly
    :param wavelength_min: Minimum wavelength that can be evaluated, in microns.
    :type wavelength_min: float, readonly
    :param wavelength_max: Maximum wavelength that can be evaluated, in microns.
    :type wavelength_max: float, readonly
    :param coefs: an array-like list of coefficients for the particular formula used.
    :type coefs: list[float], readonly
    :param absorption: The (constant) absorption coefficient used to evaluate the complex refractive index.
                        Given in cm^-1.
    :type absorption: float, readonly
    """
    cdef:
        readonly double wavelength_min
        readonly double wavelength_max
        readonly double[:] coefs
        readonly double absorption #in cm^-1
        readonly int formula_id
        dispersion_curve curve
        
    def __init__(self, int formula_id, double[:] coefs, 
                  double absorption=0.0,
                  double wavelength_min=0.1, 
                  double wavelength_max=1000000.0):
        """
        Object constructor
        
        :param int formula_id: The formula identifier {0,1,2,3,5}
        """
        self.coefs = coefs
        if formula_id==1:
            self.curve = &sellmeier_1
        elif formula_id==2:
            self.curve = &sellmeier_2
        elif formula_id==3:
            self.curve = &sellmeier_3
        elif formula_id==5:
            self.curve = &sellmeier_5
        elif formula_id==0:
            self.curve = &nondispersive_0
        else:
            raise ValueError("Unknown formula id (%d)"%(formula_id,))
        self.formula_id = formula_id
        self.absorption = absorption
        self.wavelength_min = wavelength_min
        self.wavelength_max = wavelength_max
        
    cdef np_.npy_complex128[:] c_evaluate_n(self, double[:] wavelen):
        cdef:
            dispersion_curve curve=self.curve
            double[:] coefs=self.coefs
            int i
            np_.npy_complex128[:] n_out
            double n_imag #imaginary part of n due to absorption
            double wvl
        
        n_imag = 0.0001 * self.absorption / (4* M_PI)
        
        n_out = np.empty(wavelen.shape[0], dtype=np.complex128)
        for i in range(wavelen.shape[0]):
            wvl = wavelen[i] #in microns
            if (wvl < self.wavelength_min) or (wvl > self.wavelength_max):
                msg = "Wavelength (%f) outside range of dispersion curve (%f -> %f)"%(wvl, 
                                                                                      self.wavelength_min,
                                                                                      self.wavelength_max)
                raise ValueError(msg) 
            n_out[i].real = curve(wvl, coefs)
            n_out[i].imag = n_imag * wvl
            
        return n_out
    
    def evaluate_n(self, wavelen):
        """
        Calculates the complex refractive index for the given wavelengths.
        
        :param wavelen: An array-like collection of wavelength, in microns.
        :type wavelen: double[:]
        :return: The refractive index for the given wavelength.
        :rtype: complex128[:]
        """
        return np.asarray(self.c_evaluate_n(np.asarray(wavelen)))


cdef class OpaqueMaterial(InterfaceMaterial):
    """A perfect absorber i.e. it generates no rays
    """
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        pass
    
    
cdef class TransparentMaterial(InterfaceMaterial):
    """A perfect transmitter i.e. it generates an
    outgoing ray with identical direction, polarisation etc.
    to the incoming ray. It does project the polarisation
    vectors to it's S- and P-directions, however.
    """
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        cdef:
            vector_t normal
            ray_t sp_ray
        
        normal = norm_(orient.normal)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.length = INF
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id &= ~REFL_RAY
        new_rays.add_ray_c(sp_ray)


cdef class PECMaterial(InterfaceMaterial):
    """Simulates a Perfect Electrical Conductor. I.e. incident rays are reflected with 
    100% reflectivity.
    """
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, normal
            ray_t sp_ray
            complex_t cpx
            double cosTheta
        
        normal = norm_(orient.normal)
        sp_ray = convert_to_sp(in_ray[0], normal)
        cosTheta = dotprod_(normal, in_ray.direction)
        cosThetaNormal = multvs_(normal, cosTheta)
        reflected = subvv_(in_ray.direction, multvs_(cosThetaNormal, 2))
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.direction = reflected
        sp_ray.length = INF
        sp_ray.E1_amp.real = -sp_ray.E1_amp.real
        sp_ray.E1_amp.imag = -sp_ray.E1_amp.imag
        #sp_ray.E2_amp.real = -sp_ray.E2_amp.real
        #sp_ray.E2_amp.imag = -sp_ray.E2_amp.imag
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id |= REFL_RAY
        new_rays.add_ray_c(sp_ray)
        
        
cdef class PartiallyReflectiveMaterial(InterfaceMaterial):
    """
    A simple materials with a fixed reflectivity.
    
    :param reflectivity: The material power-reflectivity given as a value from 0.0. to 1.0 
    :type reflectivity: double
    """
    cdef:
        ### power reflectivity, from 0. to 1.
        public double _reflectivity
        
    def __cinit__(self, **kwds):
        self.reflectivity = kwds.get("reflectivity", 0.5)
        
    property reflectivity:
        def __get__(self):
            return self._reflectivity
        
        def __set__(self, val):
            if val<0.0 or val>1.0:
                raise ValueError(f"Reflectivity must be in range 0.0 to 1.0 ({val} given).")
            self._reflectivity = val
        
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, normal
            vector_t tangent, tg2, in_direction
            ray_t sp_ray, sp_ray2
            double cosTheta, R, T
            complex_t P
            
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = sp_ray2 = convert_to_sp(in_ray[0], normal)
        cosTheta = dotprod_(normal, in_direction)            
        cosThetaNormal = multvs_(normal, cosTheta)
        
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        sp_ray2.accumulated_path = sp_ray.accumulated_path
        
        R = sqrt(self._reflectivity)
        T = sqrt(1-self._reflectivity)
        
        #Reflected ray
        reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.direction = reflected
        sp_ray.length = INF
        sp_ray.E1_amp *= R
        sp_ray.E2_amp *= R
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id |= REFL_RAY
        new_rays.add_ray_c(sp_ray)
            
        #Transmitted ray
        sp_ray2.origin = point
        sp_ray2.normal = normal
        sp_ray2.direction = in_direction
        sp_ray2.length = INF
        sp_ray2.E1_amp *= T
        sp_ray2.E2_amp *= T
        sp_ray2.parent_idx = idx
        sp_ray2.ray_type_id &= ~REFL_RAY
        new_rays.add_ray_c(sp_ray2)
        
        
cdef class LinearPolarisingMaterial(InterfaceMaterial):
    """Simulates a perfect polarising beam splitter. P-polarisation
    is 100% transmitted while S- is reflected"""
    
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, normal
            vector_t tangent, tg2, in_direction
            ray_t sp_ray, sp_ray2
            double cosTheta
            complex_t P
            
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = sp_ray2 = convert_to_sp(in_ray[0], normal)
        cosTheta = dotprod_(normal, in_direction)            
        cosThetaNormal = multvs_(normal, cosTheta)
        
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        sp_ray2.accumulated_path = sp_ray.accumulated_path
        
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
        sp_ray.ray_type_id |= REFL_RAY
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
        sp_ray2.ray_type_id &= ~REFL_RAY
        new_rays.add_ray_c(sp_ray2)
        
        
cdef class WaveplateMaterial(InterfaceMaterial):
    """An idealised optical retarder.
    
    :param retardance: The optical retardance, given in terms of numbers-of-wavelengths.
    :type retardance: double
    :param fast_axis: A vector giving the "fast" polarisation axis
    :type fast_axis: (double, double, double)
    """
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
        """Applies the retardance to the given Ray object.
        
        :param r: an input Ray object.
        :type r: raypier.core.ctracer.Ray
        :return: a new Ray instance
        :rtype: raypier.core.ctracer.Ray
        """
        cdef Ray out
        out = Ray()
        out.ray = self.apply_retardance_c(r.ray)
        return out
            
    
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t normal
            ray_t out_ray
            complex_t E1, retard
            
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        ###The "P" output of convert_to_sp with be aligned with the fast_axis
        ###The "S" output will thus be orthogonal to the fast axis
        out_ray = convert_to_sp(in_ray[0], self.fast_axis_)
        out_ray.accumulated_path += out_ray.length * out_ray.refractive_index.real
        out_ray.origin = point
        out_ray.normal = normal
        out_ray.direction = in_direction
        out_ray.length = INF
        out_ray.parent_idx = idx
        out_ray.ray_type_id &= ~REFL_RAY
        out_ray = self.apply_retardance_c(out_ray)
        
        new_rays.add_ray_c(out_ray)
    
    
cdef class DielectricMaterial(InterfaceMaterial):
    """Simulates Fresnel reflection and refraction at a
    normal dielectric interface.
    
    The surface normal is assumed to be pointing "out" of the material.
    
    :param complex n_inside: The refractive index on the inside of the material interface
    :param complex n_outside: The refractive index on the outside of the material interface
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
            
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2, in_direction, normal
            ray_t sp_ray
            complex_t cpx
            double cosTheta, n1, n2, N2, cos1
            double N2cosTheta, N2_sin2, tan_mag_sq, c2
            double cos2, Two_n1_cos1, aspect, T_p, T_s
            int flip
            
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
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
            sp_ray.ray_type_id |= REFL_RAY
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
            sp_ray.ray_type_id &= ~REFL_RAY
        new_rays.add_ray_c(sp_ray)
        
        
    cdef para_t eval_parabasal_ray_c(self,
                                     ray_t *base_ray, 
                                     vector_t direction, #incoming ray direction
                                   vector_t point, #position of intercept
                                   orientation_t orient,
                                   unsigned int ray_type_id, #bool, if True, it's a reflected ray
                                   ):
        cdef:
            vector_t cosThetaNormal, normal, out_dir
            para_t para_out
            double cosTheta, n1, n2
            int flip
        
        normal = norm_(orient.normal)
        direction = norm_(direction)
        cosTheta = dotprod_(normal, direction)
        cosThetaNormal = multvs_(normal, cosTheta)
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside_.real
            n2 = self.n_inside_.real
            flip = 1
        else:
            n1 = self.n_inside_.real
            n2 = self.n_outside_.real
            flip = -1
        
        if ray_type_id & REFL_RAY:
            out_dir = subvv_(direction, multvs_(cosThetaNormal, 2))
        else:
            tangent = subvv_(direction, cosThetaNormal)
            tg2 = multvs_(tangent, n1/n2)
            tan_mag_sq = mag_sq_(tg2)
            c2 = sqrt(1-tan_mag_sq)
            out_dir = subvv_(tg2, multvs_(normal, c2*flip))
        
        para_out.direction = out_dir 
        para_out.origin = point
        para_out.normal = normal
        para_out.length = INF
        return para_out
        

cdef class FullDielectricMaterial(DielectricMaterial):
    """Model for dielectric using full Fresnel equations 
    to give true phase and amplitude response
    
    :param double reflection_threshold: Sets the amplitude threshold for generating a reflected ray.
    :param double transmission_threshold: Sets the amplitude threshold for generating a transmitted ray.
    """
    cdef:
        public double reflection_threshold, transmission_threshold
        public double thickness #in microns
        complex_t n_coating_
        
    property n_coating:
        def __get__(self):
            return complex(self.n_coating_.real, self.n_coating_.imag)
        
        def __set__(self, v):
            v = complex(v)
            self.n_coating_.real = v.real
            self.n_coating_.imag = v.imag
    
    def __cinit__(self, **kwds):
        self.n_coating = kwds.get("n_coating", 1.0)
        self.thickness = kwds.get("thickness", 0.1)
        self.reflection_threshold = kwds.get('reflection_threshold', 0.1)
        self.transmission_threshold = kwds.get('transmission_threshold', 0.1)
    
    @cython.cdivision(True)
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2, in_direction, normal
            ray_t sp_ray
            complex_t cpx
            double complex n1, n2, cos2, sin2, R_p, R_s, T_p, T_s, E1_amp, E2_amp
            double cosTheta, cos1, sin1, P_in
            double tan_mag_sq, c2
            double Two_n1_cos1, aspect
            int flip
            
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        #print "in fields:", in_ray[0].E1_amp.real, in_ray[0].E1_amp.imag, in_ray[0].E2_amp.real, in_ray[0].E2_amp.imag
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
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
        #print "P_in:", P_in, n1.real, E1_amp.real, E1_amp.imag, E2_amp.real, E2_amp.imag
        if P_in==0.0:
            return
        #print "P Incoming:", P_in
        #Fresnel equations for reflection
        R_p = -(n2*cos1 - n1*cos2)/(n2*cos1 + n1*cos2)
        R_s = -(n2*cos2 - n1*cos1)/(n2*cos2 + n1*cos1)
        
        #modify in place to get reflected amplitudes
        R_s *= E1_amp
        R_p *= E2_amp
        
        #print "P_R:", (cabs(R_s)**2 + cabs(R_p)**2)
        
        if ( n1.real*(cabs(R_s)**2 + cabs(R_p)**2)/P_in ) > self.reflection_threshold:
            reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INF
            sp_ray.E1_amp.real = R_s.real
            sp_ray.E1_amp.imag = R_s.imag
            sp_ray.E2_amp.real = -R_p.real
            sp_ray.E2_amp.imag = -R_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n1.real
            sp_ray.refractive_index.imag = n1.imag
            sp_ray.ray_type_id |= REFL_RAY
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
            sp_ray.ray_type_id &= ~REFL_RAY
            new_rays.add_ray_c(sp_ray)
            
            
cdef class SingleLayerCoatedMaterial(FullDielectricMaterial):
    """
    A material with a single-layer dielectric coating at the interface.
    
    :param complex n_coating: The complex refractive index for the coating.
    :param double thickness: The thickness of the coating in microns
    """
            
    @cython.cdivision(True)
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2, in_direction, normal
            ray_t sp_ray
            complex_t cpx
            double complex n1, n2, n3, cos2, sin2, cos3, sin3, R_p, R_s, T_p, T_s, E1_amp, E2_amp
            double complex n1cos1, n2cos2, n3cos3, n1cos2, n2cos1, n2cos3, n3cos2, ep1, ep2
            double complex M00, M01, M10, M11, phi
            double cosTheta, cos1, P_in
            double tan_mag_sq, c2
            double Two_n1_cos1, aspect
            double wavelength
            int flip
            
        wavelength = self._wavelengths[in_ray.wavelength_idx]
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        E1_amp = sp_ray.E1_amp.real + 1.0j*sp_ray.E1_amp.imag
        E2_amp = sp_ray.E2_amp.real + 1.0j*sp_ray.E2_amp.imag
        cosTheta = dotprod_(normal, in_direction)
        cos1 = fabs(cosTheta)
        sin1 = sqrt(1 - cos1*cos1)
        
        n2 = self.n_coating_.real + 1.0j*self.n_coating_.imag
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside_.real + 1.0j*self.n_outside_.imag
            n3 = self.n_inside_.real + 1.0j*self.n_inside_.imag
            flip = 1
            #print "out to in", n1, n2
        else:
            n1 = self.n_inside_.real + 1.0j*self.n_inside_.imag
            n3 = self.n_outside_.real + 1.0j*self.n_outside_.imag
            flip = -1
            #print "in to out", n1, n2
            
        #apply Snell's law. These become complex.
        sin2 = (n1*sin1/n2)
        cos2 = csqrt(1 - sin2*sin2)
        
        sin3 = (n1*sin1/n3)
        cos3 = csqrt(1 - sin3*sin3)
        
        cosThetaNormal = multvs_(normal, cosTheta)
        
        #incoming power
        P_in =  n1.real*(E1_amp.real**2 + E1_amp.imag**2 + \
                         E2_amp.real**2 + E2_amp.imag**2)
        
        if P_in==0.0:
            return

        #Fresnel equations for reflection and transmission
        n1cos1 = n1*cos1
        n2cos2 = n2*cos2
        n3cos3 = n3*cos3
        dwc = 2*M_PI*self.thickness/wavelength
        phi = -I*dwc*(n2 - sin2*sin2)/cos2
        ep1 = cexp(phi)/(4*n2cos2*n3cos3)
        ep2 = cexp(-2*phi)
        M00 = -ep1*( (n1cos1-n2cos2)*(n2cos2+n3cos3) + 
                     (n1cos1+n2cos2)*(n2cos2-n3cos3)*ep2 )
                     
        M01 = ep1*( (n1cos1-n2cos2)*(n2cos2-n3cos3)*ep2 + 
                     (n1cos1+n2cos2)*(n2cos2+n3cos3) )
        
        M10 = ep1*( (n1cos1-n2cos2)*(n2cos2-n3cos3) + 
                     (n1cos1+n2cos2)*(n2cos2+n3cos3)*ep2 )
        
        M11 = -ep1*( (n1cos1-n2cos2)*(n2cos2+n3cos3)*ep2 + 
                     (n1cos1+n2cos2)*(n2cos2-n3cos3) )
        
        R_s = -M00/M01
        T_s = M10 + M11*R_s
        
        n1cos2 = n1*cos2
        n2cos1 = n2*cos1
        n2cos3 = n2*cos3
        n3cos2 = n3*cos2
        M00 = -ep1*( (n1cos2-n2cos1)*(n2cos3+n3cos2) +
                      (n1cos2+n2cos1)*(n2cos3-n3cos2)*ep2 )
        
        M01 = ep1*( (n1cos2-n2cos1)*(n2cos3-n3cos2)*ep2 +
                      (n1cos2+n2cos1)*(n2cos3+n3cos2) )
        
        M10 = ep1*( (n1cos2-n2cos1)*(n2cos3-n3cos2) +
                      (n1cos2+n2cos1)*(n2cos3+n3cos2)*ep2 )
        
        M11 = -ep1*( (n1cos2-n2cos1)*(n2cos3+n3cos2)*ep2 +
                      (n1cos2+n2cos1)*(n2cos3-n3cos2) )
        
        R_p = -M00/M01
        T_p = M10 + M11*R_p
        
        #modify in place to get reflected amplitudes
        R_s *= E1_amp
        R_p *= E2_amp
        
        
        if ( n1.real*(cabs(R_s)**2 + cabs(R_p)**2)/P_in ) > self.reflection_threshold:
            reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INF
            sp_ray.E1_amp.real = R_s.real
            sp_ray.E1_amp.imag = R_s.imag
            sp_ray.E2_amp.real = -R_p.real
            sp_ray.E2_amp.imag = -R_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n1.real
            sp_ray.refractive_index.imag = n1.imag
            sp_ray.ray_type_id |= REFL_RAY
            new_rays.add_ray_c(sp_ray)
            
        #normal transmission            
        tangent = subvv_(in_direction, cosThetaNormal)
        tg2 = multvs_(tangent, n1.real/n3.real) #This is an approximation for complex N
        tan_mag_sq = mag_sq_(tg2)
        c2 = sqrt(1-tan_mag_sq)
        transmitted = subvv_(tg2, multvs_(normal, c2*flip))
        
        aspect = sqrt(cos3.real/cos1)
        
        #modify in place to get transmitted amplitudes
        T_s *= (E1_amp*aspect)
        T_p *= (E2_amp*aspect)
        
        if ( n3.real*(cabs(T_s)**2 + cabs(T_p)**2)/P_in ) > self.transmission_threshold:
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = transmitted
            sp_ray.length = INF
            sp_ray.E1_amp.real = T_s.real
            sp_ray.E1_amp.imag = T_s.imag
            sp_ray.E2_amp.real = T_p.real
            sp_ray.E2_amp.imag = T_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n3.real
            sp_ray.refractive_index.imag = n3.imag
            sp_ray.ray_type_id &= ~REFL_RAY
            new_rays.add_ray_c(sp_ray)
    
vacuum = BaseDispersionCurve(0,np.array([1.0,]))
    
            
cdef class CoatedDispersiveMaterial(InterfaceMaterial):
    """
    A full implementation of the Fresnel Equations for a single-layer 
    coated dielectric interface. The refractive index for each ray is obtained
    by look up of the provided dispersio-curve objects for the materials on
    each side of the interface.
    
    :param dispersion_inside: The dispersion curve for the inside side of the interface
    :type dispersion_inside: raypier.core.cmaterials.BaseDispersionCurve
    :param dispersion_outside: The dispersion curve for the "outside" of the interface
    :type dispersion_outside: raypier.core.cmaterials.BaseDispersionCurve
    :param dispersion_coating: The dispersion curve for the coating material
    :type dispersion_coating: raypier.core.cmaterials.BaseDispersionCurve
    :param double coating_thickness: The coating thickness, in microns
    :param double reflection_threshold: Sets the amplitude threshold for generating a reflected ray.
    :param double transmission_threshold: Sets the amplitude threshold for generating a transmitted ray.
    """
    cdef:
        public np_.npy_complex128[:] n_inside, n_outside, n_coating
        public BaseDispersionCurve dispersion_inside
        public BaseDispersionCurve dispersion_outside
        public BaseDispersionCurve dispersion_coating
        public double coating_thickness #in microns
        public double reflection_threshold, transmission_threshold
        
    def __cinit__(self, **kwds):
        self.reflection_threshold = kwds.get('reflection_threshold', 0.1)
        self.transmission_threshold = kwds.get('transmission_threshold', 0.1)
        self.dispersion_inside = kwds.get("dispersion_inside",vacuum)
        self.dispersion_outside = kwds.get("dispersion_outside",vacuum)
        self.dispersion_coating = kwds.get("dispersion_coating",vacuum)
        self.coating_thickness = kwds.get("coating_thickness", 0.1)
        
    cdef on_set_wavelengths(self):
        cdef:
            double[:] wavelengths = self._wavelengths
        self.n_inside = self.dispersion_inside.c_evaluate_n(wavelengths)
        self.n_outside = self.dispersion_outside.c_evaluate_n(wavelengths)
        self.n_coating = self.dispersion_coating.c_evaluate_n(wavelengths)
    
    @cython.cdivision(True)
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t cosThetaNormal, reflected, transmitted
            vector_t tangent, tg2, in_direction, normal
            ray_t sp_ray
            complex_t cpx
            double complex n1, n2, n3, cos2, sin2, cos3, sin3, R_p, R_s, T_p, T_s, E1_amp, E2_amp
            double complex n1cos1, n2cos2, n3cos3, n1cos2, n2cos1, n2cos3, n3cos2, ep1, ep2
            double complex M00, M01, M10, M11, phi
            double cosTheta, cos1, P_in
            double tan_mag_sq, c2
            double Two_n1_cos1, aspect
            double wavelength
            int flip
            np_.npy_complex128 ctemp
            
        wavelength = self._wavelengths[in_ray.wavelength_idx]
        normal = norm_(orient.normal)
        in_direction = norm_(in_ray.direction)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        
        E1_amp = sp_ray.E1_amp.real + 1.0j*sp_ray.E1_amp.imag
        E2_amp = sp_ray.E2_amp.real + 1.0j*sp_ray.E2_amp.imag
        cosTheta = dotprod_(normal, in_direction)
        cos1 = fabs(cosTheta)
        sin1 = sqrt(1 - cos1*cos1)
        
        ctemp = self.n_coating[in_ray.wavelength_idx]
        n2 = ctemp.real + I*ctemp.imag#self.n_coating_.real + 1.0j*self.n_coating_.imag

        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            ctemp = self.n_outside[in_ray.wavelength_idx]
            n1 = ctemp.real + I*ctemp.imag#self.n_outside_.real + 1.0j*self.n_outside_.imag
            ctemp = self.n_inside[in_ray.wavelength_idx]
            n3 = ctemp.real + I*ctemp.imag#self.n_inside_.real + 1.0j*self.n_inside_.imag
            flip = 1
            #print( "out to in", n1, n3)
        else:
            ctemp = self.n_inside[in_ray.wavelength_idx]
            n1 = ctemp.real + I*ctemp.imag
            ctemp = self.n_outside[in_ray.wavelength_idx]
            n3 = ctemp.real + I*ctemp.imag
            flip = -1
            #print( "in to out", n1, n3 )
            
        #apply Snell's law. These become complex.
        sin2 = (n1*sin1/n2)
        cos2 = csqrt(1 - sin2*sin2)
        
        sin3 = (n1*sin1/n3)
        cos3 = csqrt(1 - sin3*sin3)
        
        cosThetaNormal = multvs_(normal, cosTheta)
        
        #incoming power
        P_in =  n1.real*(E1_amp.real**2 + E1_amp.imag**2 + \
                         E2_amp.real**2 + E2_amp.imag**2)
        
        if P_in==0.0:
            return

        #Fresnel equations for reflection and transmission
        n1cos1 = n1*cos1
        n2cos2 = n2*cos2
        n3cos3 = n3*cos3
        dwc = 2*M_PI*self.coating_thickness/wavelength
        phi = -I*dwc*(n2 - sin2*sin2)/cos2
        ep1 = cexp(phi)/(4*n2cos2*n3cos3)
        ep2 = cexp(-2*phi)
        M00 = -ep1*( (n1cos1-n2cos2)*(n2cos2+n3cos3) + 
                     (n1cos1+n2cos2)*(n2cos2-n3cos3)*ep2 )
                     
        M01 = ep1*( (n1cos1-n2cos2)*(n2cos2-n3cos3)*ep2 + 
                     (n1cos1+n2cos2)*(n2cos2+n3cos3) )
        
        M10 = ep1*( (n1cos1-n2cos2)*(n2cos2-n3cos3) + 
                     (n1cos1+n2cos2)*(n2cos2+n3cos3)*ep2 )
        
        M11 = -ep1*( (n1cos1-n2cos2)*(n2cos2+n3cos3)*ep2 + 
                     (n1cos1+n2cos2)*(n2cos2-n3cos3) )
        
        R_s = -M00/M01
        T_s = M10 + M11*R_s
        
        n1cos2 = n1*cos2
        n2cos1 = n2*cos1
        n2cos3 = n2*cos3
        n3cos2 = n3*cos2
        M00 = -ep1*( (n1cos2-n2cos1)*(n2cos3+n3cos2) +
                      (n1cos2+n2cos1)*(n2cos3-n3cos2)*ep2 )
        
        M01 = ep1*( (n1cos2-n2cos1)*(n2cos3-n3cos2)*ep2 +
                      (n1cos2+n2cos1)*(n2cos3+n3cos2) )
        
        M10 = ep1*( (n1cos2-n2cos1)*(n2cos3-n3cos2) +
                      (n1cos2+n2cos1)*(n2cos3+n3cos2)*ep2 )
        
        M11 = -ep1*( (n1cos2-n2cos1)*(n2cos3+n3cos2)*ep2 +
                      (n1cos2+n2cos1)*(n2cos3-n3cos2) )
        
        R_p = -M00/M01
        T_p = M10 + M11*R_p
        
        #modify in place to get reflected amplitudes
        R_s *= E1_amp
        R_p *= E2_amp
        
        
        if ( n1.real*(cabs(R_s)**2 + cabs(R_p)**2)/P_in ) > self.reflection_threshold:
            reflected = subvv_(in_direction, multvs_(cosThetaNormal, 2))
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = reflected
            sp_ray.length = INF
            sp_ray.E1_amp.real = R_s.real
            sp_ray.E1_amp.imag = R_s.imag
            sp_ray.E2_amp.real = -R_p.real
            sp_ray.E2_amp.imag = -R_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n1.real
            sp_ray.refractive_index.imag = n1.imag
            sp_ray.ray_type_id |= REFL_RAY
            new_rays.add_ray_c(sp_ray)
            
        #normal transmission            
        tangent = subvv_(in_direction, cosThetaNormal)
        tg2 = multvs_(tangent, n1.real/n3.real) #This is an approximation for complex N
        tan_mag_sq = mag_sq_(tg2)
        c2 = sqrt(1-tan_mag_sq)
        transmitted = subvv_(tg2, multvs_(normal, c2*flip))
        
        aspect = sqrt(cos3.real/cos1)
        
        #modify in place to get transmitted amplitudes
        T_s *= (E1_amp*aspect)
        T_p *= (E2_amp*aspect)
        
        if ( n3.real*(cabs(T_s)**2 + cabs(T_p)**2)/P_in ) > self.transmission_threshold:
            sp_ray.origin = point
            sp_ray.normal = normal
            sp_ray.direction = transmitted
            sp_ray.length = INF
            sp_ray.E1_amp.real = T_s.real
            sp_ray.E1_amp.imag = T_s.imag
            sp_ray.E2_amp.real = T_p.real
            sp_ray.E2_amp.imag = T_p.imag
            sp_ray.parent_idx = idx
            sp_ray.refractive_index.real = n3.real
            sp_ray.refractive_index.imag = n3.imag
            sp_ray.ray_type_id &= ~REFL_RAY
            new_rays.add_ray_c(sp_ray)
            
    cdef para_t eval_parabasal_ray_c(self,
                                     ray_t *base_ray, 
                                     vector_t direction, #incoming ray direction
                                   vector_t point, #position of intercept
                                   orientation_t orient,
                                   unsigned int ray_type_id, #bool, if True, it's a reflected ray
                                   ):
        cdef:
            vector_t cosThetaNormal, normal, out_dir
            para_t para_out
            double cosTheta, n1, n2
            int flip
        
        normal = norm_(orient.normal)
        direction = norm_(direction)
        cosTheta = dotprod_(normal, direction)
        cosThetaNormal = multvs_(normal, cosTheta)
        
        if cosTheta < 0.0: 
            #ray incident from outside going inwards
            n1 = self.n_outside[base_ray.wavelength_idx].real#self.n_outside_.real
            n2 = self.n_inside[base_ray.wavelength_idx].real
            flip = 1
        else:
            n1 = self.n_inside[base_ray.wavelength_idx].real
            n2 = self.n_outside[base_ray.wavelength_idx].real
            flip = -1
        
        if ray_type_id & REFL_RAY:
            out_dir = subvv_(direction, multvs_(cosThetaNormal, 2))
        else:
            tangent = subvv_(direction, cosThetaNormal)
            tg2 = multvs_(tangent, n1/n2)
            tan_mag_sq = mag_sq_(tg2)
            c2 = sqrt(1-tan_mag_sq)
            out_dir = subvv_(tg2, multvs_(normal, c2*flip))
        
        para_out.direction = out_dir 
        para_out.origin = point
        para_out.normal = normal
        para_out.length = INF
        return para_out
            

cdef class DiffractionGratingMaterial(InterfaceMaterial):
    """
    A specialised material modelling the behaviour of a reflective diffraction grating.
    
    :param double lines_per_mm: The grating line density
    :param int order: The order-of-diffraction
    :param double efficiency: A value from 0.0 to 1.0 giving the reflection efficiency of the grating
    :param origin: A vector (3-tuple) indicating the origin of the grating. Unimportant in most cases. 
                    This affects the "grating phase" which is unobservable in most situations.
    :type origin: (double, double, double)
    
    """
    ###A perfect reflection diffraction grating
    cdef:
        public double lines_per_mm #: The line-density of the grating
        public int order #order of diffraction
        public double efficiency
        vector_t origin_
        
    def __cinit__(self, **kwds):
        self.lines_per_mm = kwds.get("lines_per_mm", 1000)
        self.order = kwds.get("order", 1)
        self.efficiency = kwds.get("efficiency", 1.0)
        self.origin = kwds.get("origin", (0.0,0.0,0.0))
        
    property origin: 
        def __get__(self):
            cdef vector_t o = self.origin_
            return (o.x, o.y, o.z)
        
        def __set__(self, o):
            self.origin_.x = o[0]
            self.origin_.y = o[1]
            self.origin_.z = o[2]
    
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
           ray - the ingoing ray
           idx - the index of ray in it's RayCollection
           point - the position of the intersection (in global coords)
           normal - the outward normal vector for the surface
        """
        cdef:
            vector_t reflected, normal, tangent, tangent2
            ray_t sp_ray
            complex_t n_ray
            double k_x, k_y, k_z
            double wavelen, line_spacing
            int sign
        
        #Create surface unit axes
        normal = norm_(orient.normal)
        tangent = norm_(orient.tangent) #this is the grating wavevector direction
        tangent2 = cross_(normal, tangent) #parallel to grating lines
        
        wavelen = self._wavelengths[in_ray.wavelength_idx]
        line_spacing = 1000.0 / self.lines_per_mm #in microns
        
        reflected = norm_(in_ray.direction) #ensure unit vector
        #Break incoming direction vector into components along surface unit axes
        k_z = dotprod_(normal, reflected)
        k_y = dotprod_(tangent2, reflected)
        k_x = dotprod_(tangent, reflected)
        
        if k_z<0.0:#invert for a reflecting grating
            sign = 1
        else:
            sign = -1
        
        n_ray = in_ray[0].refractive_index
        k_x = k_x - self.order*wavelen/(line_spacing*n_ray.real)
        #y-component (parallel to grating lines) doesn't change
        
        k_z = 1 - (k_x*k_x) - (k_y*k_y) #Apply Pythagoras to get z-component
        if k_z < 0: #diffracted ray is evanescent
            return
        k_z = sign * sqrt(k_z) 
        
        #Add the x- , y- and z-components together
        reflected = multvs_(tangent, k_x)
        reflected = addvv_(reflected, multvs_(tangent2, k_y))
        reflected = addvv_(reflected, multvs_(normal, k_z))
                
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real
        
        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.direction = reflected
        sp_ray.E1_amp.real = -sp_ray.E1_amp.real * self.efficiency
        sp_ray.E1_amp.imag = -sp_ray.E1_amp.imag * self.efficiency
        sp_ray.E2_amp.real = sp_ray.E2_amp.real * self.efficiency
        sp_ray.E2_amp.imag = sp_ray.E2_amp.imag * self.efficiency
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id |= REFL_RAY
        
        ### This is the mysterious Grating Phase
        #Need to convert origin offset to microns since line-spacing is in microns
        sp_ray.phase += 1000.0*dotprod_(subvv_(self.origin_, point), tangent)*self.order*2*M_PI/line_spacing
        
        new_rays.add_ray_c(sp_ray)
        
    cdef para_t eval_parabasal_ray_c(self,
                                     ray_t *base_ray, 
                                     vector_t direction, #incoming ray direction
                                   vector_t point, #position of intercept
                                   orientation_t orient,
                                   unsigned int ray_type_id, #bool, if True, it's a reflected ray
                                   ):        
        cdef:
            para_t para_out
            vector_t reflected, normal, tangent, tangent2
            ray_t sp_ray
            complex_t n_ray
            double k_x, k_y, k_z
            double wavelen, line_spacing
            int sign
        
        #Create surface unit axes
        normal = norm_(orient.normal)
        tangent = norm_(orient.tangent) #this is the grating wavevector direction
        tangent2 = cross_(normal, tangent) #parallel to grating lines
        
        wavelen = self._wavelengths[base_ray.wavelength_idx]
        line_spacing = 1000.0 / self.lines_per_mm #in microns
        
        reflected = norm_(direction) #ensure unit vector
        #Break incoming direction vector into components along surface unit axes
        k_z = dotprod_(normal, reflected)
        k_y = dotprod_(tangent2, reflected)
        k_x = dotprod_(tangent, reflected)
        
        if k_z<0.0:#invert for a reflecting grating
            sign = 1
        else:
            sign = -1
        
        n_ray = base_ray.refractive_index #Should this be the parent refractive index?
        k_x = k_x - self.order*wavelen/(line_spacing*n_ray.real)
        #y-component (parallel to grating lines) doesn't change
        
        k_z = 1 - (k_x*k_x) - (k_y*k_y) #Apply Pythagoras to get z-component
        if k_z < 0: #diffracted ray is evanescent
            print("Error. Parabasal reflection is imaginary!")
            #return #Then we're in trouble!
        k_z = sign * sqrt(k_z) 
        
        #Add the x- , y- and z-components together
        reflected = multvs_(tangent, k_x)
        reflected = addvv_(reflected, multvs_(tangent2, k_y))
        reflected = addvv_(reflected, multvs_(normal, k_z))
                
        para_out.direction = reflected
        para_out.origin = point
        para_out.normal = normal
        para_out.length=INF        
        
        return para_out 
        
        
cdef class CircularApertureMaterial(InterfaceMaterial):
    """Similar to the TransparentMaterial i.e. it generates an
    outgoing ray with identical direction, polarisation etc.
    to the incoming ray. The material attenuates the E_field amplitudes
    according to the radial distance from the surface origin. 
    
    :param double outer_radius: Rays passing outside this radius are not intercepted.
    :param double radius: The radius of the inner hole.
    :param double edge_width: Rays passing through the inner hole are attenuated according to their 
                    proximity to the edge. The edge-width sets the width of the error-function (errf)
                    used to calculate the attenuation.
    :param int invert: Inverts the aperture to make a field-stop.
    :param origin: The centre point of the aperture.
    :type origin: (double, double, double)
    """
    cdef:
        public double outer_radius
        public double radius
        public double edge_width
        public int invert
        vector_t origin_
        
    def __cinit__(self, **kwds):
        self.outer_radius = kwds.get("outer_radius", 25.0)
        self.radius = kwds.get("radius", 15.0)
        self.edge_width = kwds.get("edge_width", 1.0)
        self.origin = kwds.get("origin", (0.0,0.0,0.0))
        self.invert = kwds.get("invert", 0)
        
    property origin: 
        def __get__(self):
            cdef vector_t o = self.origin_
            return (o.x, o.y, o.z)
        
        def __set__(self, o):
            self.origin_.x = o[0]
            self.origin_.y = o[1]
            self.origin_.z = o[2]
        
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        cdef:
            vector_t normal
            ray_t sp_ray
            double atten, r
            double width = self.edge_width
            
        r = sqrt(mag_sq_(subvv_(self.origin_, point)))
        if (r > self.outer_radius):
            return
        atten = 0.5 + 0.5*erf((self.radius - r)/width)
        if self.invert:
            atten = 1 - atten
        
        normal = norm_(orient.normal)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real

        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id &= ~REFL_RAY
        
        sp_ray.E1_amp.real *= atten
        sp_ray.E1_amp.imag *= atten
        sp_ray.E2_amp.real *= atten
        sp_ray.E2_amp.imag *= atten
        
        new_rays.add_ray_c(sp_ray)


cdef class RectangularApertureMaterial(InterfaceMaterial):
    """
    A rectangular aperture.
    
    :param double outer_width:
    :param double outer_height:
    :param double width:
    :param double height:
    :param double edge_width:
    :param int invert:
    :param origin:
    :type origin: (double, double, double)
    """
    cdef:
        public double outer_width
        public double outer_height
        public double width
        public double height
        public double edge_width
        public int invert
        vector_t origin_
        
    def __cinit__(self, **kwds):
        self.outer_width = kwds.get("outer_width", 15.0)
        self.outer_height = kwds.get("outer_height", 20.0)
        self.width = kwds.get("width", 5.0)
        self.height = kwds.get("height", 10.0)
        self.edge_width = kwds.get("edge_width", 1.0)
        self.origin = kwds.get("origin", (0.0,0.0,0.0))
        self.invert = kwds.get("invert", 0)
        
    property origin: 
        def __get__(self):
            cdef vector_t o = self.origin_
            return (o.x, o.y, o.z)
        
        def __set__(self, o):
            self.origin_.x = o[0]
            self.origin_.y = o[1]
            self.origin_.z = o[2]
        
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        cdef:
            vector_t normal, p
            ray_t sp_ray
            double atten, px, py
            double width = self.edge_width
            double x = self.width/2.
            double y = self.height/2.
            
        p = subvv_(point, self.origin_)
        px = dotprod_(p,orient.tangent)
        py = dotprod_(p,cross_(orient.normal,orient.tangent))
            
        if fabs(px) > self.outer_width/2.:
            return 
        if fabs(py) > self.outer_height/2.:
            return 
            
        atten = 0.5 - 0.5*erf((px - x)/width)
        atten *= 0.5 - 0.5*erf(-(px + x)/width)
        atten *= 0.5 - 0.5*erf((py - y)/width)
        atten *= 0.5 - 0.5*erf(-(py + y)/width)
        
        if self.invert:
            atten = 1 - atten
        
        normal = norm_(orient.normal)
        sp_ray = convert_to_sp(in_ray[0], normal)
        sp_ray.accumulated_path += sp_ray.length * sp_ray.refractive_index.real

        sp_ray.origin = point
        sp_ray.normal = normal
        sp_ray.parent_idx = idx
        sp_ray.ray_type_id &= ~REFL_RAY
        
        sp_ray.E1_amp.real *= atten
        sp_ray.E1_amp.imag *= atten
        sp_ray.E2_amp.real *= atten
        sp_ray.E2_amp.imag *= atten
        
        new_rays.add_ray_c(sp_ray)
        
        
cdef class ResampleGaussletMaterial(InterfaceMaterial):
    """
    This is a special pseudo-material which generates new rays not by
    the normal process of refraction or reflection of an incoming ray,
    but by computing a new set of Gausslets by computing the E-field
    at a set of grid points and launching new Gausslets from these points.
    
    The material needs the set of new launch-positions to be given up front (i.e. before tracing). 
    
    :param int size: The size of the new GaussletCollection to allocate up front.
    :param callable eval_func: A callable taking a GaussletCollection as its single argument.
            The new outgoing rays will be returned by this callable.
    """
    cdef:
        public int capture_count
        public GaussletCollection captured_rays
        object _evaluation_func
        
    def __cinit__(self, **kwds):
        size = kwds.get("size", 2)
        func = kwds.get("eval_func", None)
        self.captured_rays = GaussletCollection(size)
        if func is not None:
            self.eval_func = func
        
    property eval_func:
        def __get__(self):
            return self._evaluation_func
        
        def __set__(self, obj):
            if not callable(obj):
                raise ValueError("The eval_func property must be a callable.")
            self._evaluation_func = obj
            
    def is_decomp_material(self):
        return True
    
    cdef void eval_child_ray_c(self,
                            ray_t *in_ray, 
                            unsigned int idx, 
                            vector_t point,
                            orientation_t orient,
                            RayCollection new_rays):
        """
        Capture each incident ray in an internal GaussletCollection.
        """
        cdef:
            gausslet_t *in_g
            
        if in_ray[0].ray_type_id & 2: #check GAUSSLET bit
            in_g = <gausslet_t *>in_ray #oooh, dangerous! Works because each gausslet_t starts with a ray_t
            self.captured_rays.add_gausslet_c(in_g[0])
            self.capture_count += 1
    
    cdef void eval_decomposed_rays_c(self, GaussletCollection new_rays):
        """
        Here we compute the new rays and append them to the new_rays object.
        """
        cdef:
            GaussletCollection gc
            
        self.captured_rays._wavelengths = new_rays._wavelengths
        gc = self._evaluation_func(self.captured_rays)
        new_rays.extend_c(gc)
        #Clear the captured rays
        self.captured_rays.n_rays = 0

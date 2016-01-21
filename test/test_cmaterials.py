#!/usr/bin/env python

import pyximport
pyximport.install()

import sys
sys.path.append('..')
from raytrace.cmaterials import FullDielectricMaterial, Convert_to_SP, \
    TransparentMaterial
from raytrace.ctracer import Ray, RayCollection, norm, dotprod, subvv, cross
import unittest
from math import sqrt
import random
import numpy


P = lambda z: z.real**2 + z.imag**2

def ray_power(ray):
    """incident power per unit area"""
    P1 = (ray.E1_amp.real**2 + ray.E1_amp.imag**2)*ray.refractive_index.real
    P2 = (ray.E2_amp.real**2 + ray.E2_amp.imag**2)*ray.refractive_index.real
    ### We don't need to handle incident aspect area here as the ray already compensate for this
    #aspect = abs(dotprod(norm(ray.direction), normal))
    return (P1+P2) 


class TestTransparentMaterial(unittest.TestCase):
    def test_preserve_polariasation(self):
        m = TransparentMaterial
        ray_idx = 123
        ray_in = Ray(origin=(0.0,0.0,0.0),
                     direction=(0.0,0.0,1.0),
                     E_vector=(1.0,0.0,0.0),
                     E1_amp = (1.0+2.0j),
                     E2_amp = (3.0+4.0j))
        
        out = RayCollection(1)
        point = (0.0,0.0,0.0)
        normal = (0.0,0.0,-1.0)
        
        
        


class TestConvertToSP(unittest.TestCase):
    def test_cross(self):
        v1 = (1.0,0.0,0.0)
        v2 = (0.0,1.0,0.0)
        self.assertEquals(cross(cross(v1,v2),v1), v2)
        
    
    def test_project_E(self):
        ray_in = Ray(origin=(-1.0,-2.0,-3.0),
                     direction=(1.0,2.0,3.0),
                     E_vector=(1.0,-2.0,0.0),
                     E1_amp = (1.0+2.0j),
                     E2_amp = (3.0+4.0j))
        ray_in.project_E(1.0,-2.0,0.0)
        
        E_vector = ray_in.E_vector
        E1_amp = ray_in.E1_amp
        E2_amp = ray_in.E2_amp
        
        ray_in.project_E(0.3,.7,0)
        ray_in.project_E(-0.8,0.235,0)
        ray_in.project_E(1.0,-2.0,0.0)
        
        self.assertEquals(E_vector, ray_in.E_vector)
        self.assertAlmostEquals(E1_amp, ray_in.E1_amp)
        self.assertAlmostEquals(E2_amp, ray_in.E2_amp)
        
        
    
    def test_conserve_power(self):
        ray_in = Ray(origin=(-1.0,-2.0,-3.0),
                     direction=(1.0,2.0,3.0),
                     E_vector=(1.0,-2.0,0.0),
                     E1_amp = (1.0+2.0j),
                     E2_amp = (3.0+4.0j))
        
        P_in = P(ray_in.E1_amp) + P(ray_in.E2_amp)
        
        normal = (0.2,0.3,-1.1)
        
        ray_out = Convert_to_SP(ray_in, normal)
        
        P_out = P(ray_out.E1_amp) + P(ray_out.E2_amp)
        
        self.assertAlmostEqual(P_in, P_out)
        
        ray_out.project_E(1.0,-2.0,0.0)
        
        self.assertAlmostEquals(ray_out.E1_amp, ray_in.E1_amp)
        
    def test_normal_incidence(self):
        ray_in = Ray(origin=(-1.0,-2.0,-3.0),
                     direction=(1.0,2.0,3.0),
                     E_vector=(1.0,-2.0,0.0),
                     E1_amp = (1.0+2.0j),
                     E2_amp = (3.0+4.0j))
        
        P_in = P(ray_in.E1_amp) + P(ray_in.E2_amp)
        
        normal = (1.0,2.0,3.0)
        
        ray_out = Convert_to_SP(ray_in, normal)
        
        P_out = P(ray_out.E1_amp) + P(ray_out.E2_amp)
        
        self.assertAlmostEqual(P_in, P_out)
        
        ray_in.project_E(1.0,-2.0,0.0)
        ray_out.project_E(1.0,-2.0,0.0)
        
        self.assertEquals(ray_out.E_vector, ray_in.E_vector)
        self.assertAlmostEquals(ray_out.E1_amp, ray_in.E1_amp)
        self.assertAlmostEquals(ray_out.E2_amp, ray_in.E2_amp)


class TestFullDielectricMaterial(unittest.TestCase):
    def setUp(self):
        pass
    
    def test_eval_child_ray(self):
        ray_in = Ray(origin=(-1.0,-2.0,-3.0),
                     direction=(1.0,2.0,3.0),
                     E_vector=(1.0,-2.0,0.0),
                     E1_amp = (1.0+1.0j),
                     E2_amp = (2.0+0.0j),
                     refractive_index=1.5)
        
        mat = FullDielectricMaterial(n_inside=1.0,
                                     n_outside=1.5,
                                     reflection_threshold=-0.01,
                                     transmission_threshold=-0.01)
        
        out = RayCollection(8)
        ray_idx = 123
        point = (0.0,0.0,0.0)
        normal = (0.0,0.0,-1.0)
        
        mat.eval_child_ray(ray_in, 
                           ray_idx,
                           point,
                           normal,
                           out)
        
        self.assertEquals(len(out),2)
        
        aspect_in = dotprod(norm(ray_in.direction), norm(normal))
        P_in = (P(ray_in.E1_amp) + P(ray_in.E2_amp))*(ray_in.refractive_index.real)
        
        print "D_in:", norm((1.0,2.0,3.0)), "P_in:", P_in*aspect_in
        for item in out:
            E1 = item.E1_amp
            E2 = item.E2_amp
            aspect_out = dotprod(norm(item.direction), norm(normal))
            print (P(E1) + P(E2))*(item.refractive_index.real)*aspect_out, item.refractive_index
        
    def test_normal_incidence(self):
        n_out = 1.3
        n_in = 2.0
        
        ray_in = Ray(origin=(0.0,0.0,-1.0),
                     direction=(0.0,0.0,3.0),
                     E_vector=(5.0,0.0,0.0),
                     E1_amp = (1.0+1.0j),
                     E2_amp = (2.0+0.0j),
                     refractive_index=n_out)
        
        mat = FullDielectricMaterial(n_inside=n_in,
                                     n_outside=n_out,
                                     reflection_threshold=-0.01,
                                     transmission_threshold=-0.01)
        
        out = RayCollection(8)
        ray_idx = 123
        point = (0.0,0.0,0.0)
        normal = (0.0,0.0,-1.0)
        
        P_in = ray_power(ray_in, normal)
        
        mat.eval_child_ray(ray_in, 
                           ray_idx,
                           point,
                           normal,
                           out)
        
        self.assertEquals(len(out),2)
        
        #textbook reflection and transmission coefficients
        R = ((n_in - n_out)/(n_in + n_out))**2
        T = (n_in/n_out)*((2*n_out/(n_out+n_in))**2)
        
        self.assertAlmostEquals(1, T+R)
    
        refl_pow = ray_power(out[0], normal)
        self.assertAlmostEquals(R, refl_pow/P_in)
        
        trans_pow = ray_power(out[1], normal)
        self.assertAlmostEquals(T, trans_pow/P_in)
        
        self.assertAlmostEqual(P_in, refl_pow+trans_pow)
        
    def test_brewster_angle(self):
        n_out = 1.2
        n_in = 2.0
        
        theta_B = numpy.arctan2(n_in, n_out)
        print "Brewsters angle:", theta_B*180/numpy.pi
        
        y = numpy.sin(theta_B)
        z = numpy.cos(theta_B)
        
        #angle inside dielectric (for transmitted ray)
        theta_inside = numpy.arcsin( n_out * numpy.sin(theta_B) / n_in )
        trans_direction = norm((0.0, numpy.sin(theta_inside), numpy.cos(theta_inside)))
        
        ###input ray is directed into object
        ray_in = Ray(origin=(0.0,-y,-z),
                     direction=(0.0,y,z),
                     E_vector=(0.0,1.0,0.0),
                     E1_amp = (1.0+1.0j),
                     E2_amp = (0.0+0.0j),
                     refractive_index=n_out)
        
        mat = FullDielectricMaterial(n_inside=n_in,
                                     n_outside=n_out,
                                     reflection_threshold=-0.01,
                                     transmission_threshold=-0.01)
        
        out = RayCollection(8)
        ray_idx = 123
        point = (0.0,0.0,0.0)
        normal = (0.0,0.0,-1.0)
        
        P_in = ray_power(ray_in)
        
        mat.eval_child_ray(ray_in, 
                           ray_idx,
                           point,
                           normal,
                           out)
        
        self.assertEquals(len(out),2)
        
        refl_pow = ray_power(out[0])
        trans_pow = ray_power(out[1])
        
        self.assertAlmostEquals(0.0, refl_pow)
        self.assertAlmostEquals(refl_pow, out[0].power)
        
        self.assertAlmostEquals(P_in, trans_pow)
        
        ###check relfected ray direction
        self.assertAlmostEquals(abs(dotprod(norm(subvv(ray_in.direction, 
                                                       out[0].direction)),
                                            normal)),
                                1)
        
        ###check transmitted ray direction
        self.assertAlmostEquals(abs(dotprod(out[1].direction,
                                            trans_direction)), 1)
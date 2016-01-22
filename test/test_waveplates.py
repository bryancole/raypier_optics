
import unittest

from raytrace.tracer import RayTraceModel
from raytrace.waveplates import Waveplate
from raytrace.sources import ParallelRaySource

from numpy import pi
import numpy


class TestWaveplate(unittest.TestCase):
    def setUp(self):
        self.wp = Waveplate(direction=(0.,0.,1.))
        self.src = ParallelRaySource(rings=0,
                                     E_vector=(1,0,0),
                                     E1_amp = 1.0+0.0j,
                                     E2_amp = 0.0+0.0j,
                                     direction = (0,0,1),
                                     origin = (0,0,-10),
                                    )
        
        self.model = RayTraceModel(optics=[self.wp,],
                                   sources=[self.src])
        
    def test_retardance_property(self):
        wp = self.wp
        for val in [0.0, 0.25, 0.5, 0.75, 1.0]:
            wp.retardance = val
            self.assertEquals(wp.retardance, val)
        
    def test_input_ray(self):
        ray = self.src.InputRays[0]
        print ray.power, ray.ellipticity, ray.E_vector, ray.E1_amp, ray.E2_amp
        self.assertEquals(ray.ellipticity, 0.0)
        
    def test_zero_retardance(self):
        wp = self.wp
        wp.retardance = 0.0
        wp.rotation = 19.235
        print wp.x_axis
        self.model.trace_all()
        
        ray_in = self.src.InputRays[0]
        ray_out = self.src.TracedRays[-1][0]

        ray_in.project_E(1,1,0)
        ray_out.project_E(1,1,0)
        
        self.assertEquals(ray_in.E_vector, ray_out.E_vector)
        self.assertAlmostEquals(ray_in.E1_amp, ray_out.E1_amp)
        self.assertAlmostEquals(ray_in.E2_amp, ray_out.E2_amp)
        
    def test_half_wave_plate(self):
        wp = self.wp
        wp.retardance = 0.5
        for rot in numpy.linspace(-180,180,20):
            wp.rotation = rot
            self.model.trace_all()
            
            ray_in = self.src.InputRays[0]
            ray_out = self.src.TracedRays[-1][0]
    
            self.assertAlmostEquals(ray_in.ellipticity, ray_out.ellipticity)
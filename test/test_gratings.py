
import unittest

import pyximport
pyximport.install()

from raypier.sources import SingleRaySource
from raypier.diffraction_gratings import RectangularGrating
from raypier.tracer import RayTraceModel
from raypier.core.ctracer import cross

import numpy


class TestGratingEquation(unittest.TestCase):
    def test_grating_equation1(self):
        grating = RectangularGrating(centre=(0,0,0),
                                     direction=(0,1,0),
                                     order=-1,
                                     lines_per_mm=1400)
        source = SingleRaySource(origin=(-10,10,0),
                                 direction=(1,-1,0),
                                 wavelength=0.8,
                                 )
        model = RayTraceModel(sources=[source,],
                              optics=[grating,])
        model.trace_all()
        traced_rays = list(source.traced_rays)
        self.assertEqual(len(traced_rays),2)
        
        out = traced_rays[-1]
        in_ray = traced_rays[0][0]
        out_ray = out[0]
        
        
        sin_in = -cross(in_ray.direction, grating.direction)[2]
        sin_out = cross(out_ray.direction, grating.direction)[2]
        
        m = grating.order
        D=1000./grating.lines_per_mm #in microns
        wl=0.8 #wavelength in microns
        LHS = m*wl-(D*sin_in)
        RHS = D*sin_out
        self.assertAlmostEqual(abs(LHS), abs(RHS))
        
    def test_grating_phase(self):
        grating = RectangularGrating(centre=(0.3e-3,0.3e-3,0.4e-3),
                                     direction=(0,1,0),
                                     order=-1,
                                     lines_per_mm=1400)
        source = SingleRaySource(origin=(-10,10,0),
                                 direction=(1,-1,0),
                                 wavelength=0.8,
                                 )
        model = RayTraceModel(sources=[source,],
                              optics=[grating,])
        model.trace_all()
        traced_rays = list(source.traced_rays)
        self.assertEqual(len(traced_rays),2)
        
        out = traced_rays[-1]
        in_ray = traced_rays[0][0]
        out_ray = out[0]
        
        print(out_ray.phase)

        

import unittest
import numpy

from raytrace.sources import ConfocalRayFieldSource


class TestConfocalPhase(unittest.TestCase):
    def test_phase_at_focus(self):
        
        src = ConfocalRayFieldSource(angle=5.0,
                                     angle_step=0.5,
                                     origin=(0,0,0),
                                     working_dist=100.0,
                                     direction=(0,0,1),
                                     wavelength=100.0
                                     )
        rays = src.InputRays.copy_as_array()
        
        vec_to_focus = rays['origin'] - numpy.array([[0,0,src.working_dist]])
        
        path_to_focus = numpy.sqrt((vec_to_focus**2).sum(axis=-1)) 
        offset = 2000.0*path_to_focus*numpy.pi/src.wavelength
        ph_at_origin = rays['phase'] + offset
        
        self.assertTrue( numpy.allclose(ph_at_origin, 0.0))
        
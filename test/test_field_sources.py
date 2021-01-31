
import unittest
import numpy

from raypier.sources import ConfocalRayFieldSource
from raypier.mirrors import PlanarWindow
from raypier.tracer import RayTraceModel
from raypier.core.ctracer import cross
from raypier.core.cmaterials import Convert_to_SP


class TestConfocalPhase(unittest.TestCase):
    def test_phase_at_focus(self):
        
        src = ConfocalRayFieldSource(angle=5.0,
                                     angle_step=0.5,
                                     origin=(0,0,0),
                                     working_dist=100.0,
                                     direction=(0,0,1),
                                     wavelength=100.0
                                     )
        rays = src.input_rays.copy_as_array()
        
        vec_to_focus = rays['origin'] - numpy.array([[0,0,src.working_dist]])
        
        path_to_focus = numpy.sqrt((vec_to_focus**2).sum(axis=-1)) 
        offset = 2000.0*path_to_focus*numpy.pi/src.wavelength
        ph_at_origin = rays['phase'] + offset
        
        self.assertTrue( numpy.allclose(ph_at_origin, 0.0))
        
    def test_EField_after_tracing(self):
        src = ConfocalRayFieldSource(angle=5.0,
                                     angle_step=4.0,
                                     origin=(0,0,0),
                                     working_dist=100.0,
                                     direction=(0,0,1),
                                     wavelength=100.0,
                                     )
        
        window = PlanarWindow(centre=(0,0,20),
                              diameter=25.0,
                              thickness=5.0,
                              n_inside=1.0)
        
        model = RayTraceModel(optics=[window],
                              sources=[src])
        model.trace_all()
        
        for rays in src.traced_rays:
            E = rays.E_vector
            H = numpy.cross(rays.direction, E)
            F = E * rays.E1_amp[:,None] + H * rays.E2_amp[:,None]
            print(F.imag)
        
        rays = src.traced_rays
        
        ray1 = rays[0][0]
        ray2 = rays[1][0]
        
        self.assertAlmostEqual(ray1.direction, ray2.direction)
        
        #print(ray2.normals)
        ray2b = Convert_to_SP(ray1, ray2.normals)
        #print(ray1.E_vector, ray2.E_vector)
        
        for ray in (ray1, ray2, ray2b):
            E = ray.E_vector
            H = cross(ray.direction, E)
            F = tuple(e*ray.E1_amp.real + h*ray.E2_amp.real for e,h in zip(E,H))
            print( F)
        
        
         
        
        
        
        
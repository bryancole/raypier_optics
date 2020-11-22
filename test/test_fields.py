
import unittest

from raytrace.sources import HexagonalRayFieldSource
from raytrace.fields import evaluate_neighbours, project_to_sphere,\
        evaluate_modes
        
from raytrace.cfields import inv_area_of_ellipse

import numpy


class TestProjectToSphere(unittest.TestCase):
    def test_project1(self):
        src = HexagonalRayFieldSource(centre=(0,0,-50),
                                      direction=(0,0,1),
                                      spacing=2,
                                      radius=10.0)
        
        rays = src.InputRays.copy_as_array()
        phase = rays['phase'] #will be zeros
        wavelengths= src.wavelength_list
        
        out_rays, out_phase = project_to_sphere(rays, phase, wavelengths, 
                                     centre=(0,0,0), 
                                     radius=20.0)
        
        print(out_rays['origin'])
        self.assertAlmostEqual(out_rays['origin'][:,2].min(), -20.0)
        
        
class TestEvalNeighbours(unittest.TestCase):
    def test_eval_neighbours(self):
        src = HexagonalRayFieldSource(centre=(0,0,-50),
                                      direction=(0,0,1),
                                      spacing=4,
                                      radius=10.0)
        
        rays_in = src.InputRays.copy_as_array()
        nb_idx = src.neighbours

        self.assertEquals(len(rays_in), 19)
        rays, x, y, dx, dy = evaluate_neighbours(rays_in, nb_idx)
        self.assertEquals(len(rays), 7)
        self.assertEqual( len(x.shape), 2)
        
        offset = x**2 + y**2
        self.assertTrue( numpy.allclose(offset,offset.mean()) )
        
        self.assertTrue( numpy.allclose(dx, 0) )
        self.assertTrue( numpy.allclose(dy, 0) )
        
        
class TestEvalModes(unittest.TestCase):
    def test_solve_modes(self):
        src = HexagonalRayFieldSource(centre=(0,0,-50),
                                      direction=(0,0,1),
                                      spacing=4,
                                      radius=10.0)
        
        rays_in = src.InputRays.copy_as_array()
        nb_idx = src.neighbours

        self.assertEquals(len(rays_in), 19)
        rays, x, y, dx, dy = evaluate_neighbours(rays_in, nb_idx)
        print(x.shape)
        Z = evaluate_modes(rays, x, y, dx, dy)
        print(Z.shape)
        
        
class TestAreaOfEllipse(unittest.TestCase):
    def test_area_of_ellipse(self):
        ra = 3.0
        rb = 8.0
        angle = 30*numpy.pi/180.0
        
        a=ra #(1/(ra**2))
        b=rb #(1/(rb**2))
        
        sin = numpy.sin(angle)
        cos = numpy.cos(angle)
        
        A = (a**2)*(sin**2) + (b**2)*(cos**2)
        B = 2*(b**2 - a**2)*sin*cos
        C = (a**2)*(cos**2) + (b**2)*(sin**2)
        
        area = numpy.pi*ra*rb
        print("True area:", area)
        
        area2 = (inv_area_of_ellipse(A, B, C)**1)
        
        self.assertAlmostEqual(area, area2)

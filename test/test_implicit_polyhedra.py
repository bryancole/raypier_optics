
import unittest
import numpy as np
print("Numpy:", np.__version__)

from raypier.implicit_polyhedra import Plane, get_intersection_of_planes, Intersection, Polyhedron



class TestPlanesIntersection(unittest.TestCase):
    def test_simple_intersection(self):
        p1 = Plane(origin=(1.2,3.2,2.7), normal=(1.0,0.,0.)) #normal is outward direction
        p2 = Plane(origin=(2.5,1.6,7.7), normal=(0.,1.,0.))
        p3 = Plane(origin=(-2.5,2.6,5.1), normal=(0.,0.,1.))
        
        print(get_intersection_of_planes(p1,p2,p3))
        
        v = p1.c_surface.evaluate(0.,0.,0.)
        print(f"eval {v}")
        v = p1.c_surface.evaluate(3.,7.,-9.)
        print(f"eval {v}")
        
        
    def test_find_vertices(self):
        p1 = Plane(origin=(1.2,3.2,2.7), normal=(-1.0,0.,0.))
        p2 = Plane(origin=(2.5,1.6,7.7), normal=(0.,-1.,0.))
        p3 = Plane(origin=(-2.5,2.6,5.1), normal=(0.,0.,-1.))
        p4 = Plane(origin=(6.5,8.6,7.1), normal=(1.,1.,1.))
        
        phdr = Polyhedron(planes=[p1,p2,p3,p4])
        
        for v in phdr.get_vertices():
            print(v)
            
    
class TestPolyhedron(unittest.TestCase):
    def test_visualisation(self):
        from raypier.tracer import RayTraceModel
        from raypier.sources import ParallelRaySource
        
        p1 = Plane(origin=(1.2,3.2,2.7), normal=(-1.0,0.,0.))
        p2 = Plane(origin=(2.5,1.6,1.7), normal=(0.,-1.,0.))
        p3 = Plane(origin=(-2.5,2.6,2.1), normal=(0.,0.,-1.))
        p4 = Plane(origin=(10.,10.,10.0), normal=(1.,1.,1.))
        
        p5 = Plane(origin=(5.,5.,5.), normal=(-1.,-1.,-1.))
        
        p6 = Plane(origin=(0.,-6.,0.), normal=(1.,-1,0))
        
        phdr = Polyhedron(planes=[p1,p2,p3,p4,p5,p6])

        src = ParallelRaySource(origin=(0,0,-20), direction=(0,0,1.))
        
        model = RayTraceModel(optics=[phdr,],sources=[src,])
        model.configure_traits()
        
    def test_view_verts(self):
        from raypier.tracer import RayTraceModel
        
        p1 = Plane(origin=(1.2,3.2,2.7), normal=(-1.0,0.,0.))
        p2 = Plane(origin=(2.5,1.6,1.7), normal=(0.,-1.,0.))
        p3 = Plane(origin=(-2.5,2.6,2.1), normal=(0.,0.,-1.))
        p4 = Plane(origin=(10.,10.,10.0), normal=(1.,1.,1.))
        p5 = Plane(origin=(5.,5.,5.), normal=(-1.,-1.,-1.))
        
        phdr = Polyhedron(planes=[p1,p2,p3,p4,p5])
        
        for vert in phdr.get_vertices():
            print(vert)
            
        
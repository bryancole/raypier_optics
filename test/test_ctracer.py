#!/usr/bin/env python

import pyximport
pyximport.install()

import sys
sys.path.append('..')
from raytrace import ctracer
import unittest
from math import sqrt


class TestVectorMathsFunctions(unittest.TestCase):
    def test_set_v(self):
        a = (1., 2., 3.)
        b = ctracer.py_set_v(a)
        self.assertEquals(a,b)
    
    def test_addvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.addvv(a,b)
        self.assertEquals(c, tuple(ai+bi for ai,bi in zip(a,b)))
        
    def test_subvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.subvv(a,b)
        self.assertEquals(c, tuple(ai-bi for ai,bi in zip(a,b)))
        
    def test_multvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.multvv(a,b)
        self.assertEquals(c, tuple(ai*bi for ai,bi in zip(a,b)))
        
    def test_addvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.addvs(a,b)
        self.assertEquals(c, tuple(ai+b for ai in a))
        
    def test_subvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.subvs(a,b)
        self.assertEquals(c, tuple(ai-b for ai in a))
        
    def test_multvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.multvs(a,b)
        self.assertEquals(c, tuple(ai*b for ai in a))
        
    def test_sep(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.sep(a,b)
        self.assertEquals(c, sqrt(sum((ai-bi)**2 for ai,bi in zip(a,b))))
        
    def test_dotprod(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.dotprod(a,b)
        self.assertEquals(c, sum(x*y for x,y in zip(a,b)))
        
    def test_cross(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.cross(a,b)
        A = ctracer.dotprod(a,c)
        B = ctracer.dotprod(b,c)
        self.assertEquals(A, 0.0)
        self.assertEquals(B, 0.0)
        
    def test_map(self):
        a = [3.4, 5.6, 7.8]
        b = ctracer.mag(a)
        self.assertEquals(b, sqrt(sum(x**2 for x in a)))
        
    def test_norm(self):
        a = [1.1, 2.2, 3.3]
        b = ctracer.norm(a)
        print "B", b
        self.assertAlmostEqual(1.0, sum(x**2 for x in b))
        self.assertAlmostEqual(0.0, ctracer.mag(ctracer.cross(a,b)))
        
        
class TestPECMaterial(unittest.TestCase):
    def test_eval_child(self):
        m = ctracer.PECMaterial()
        in_ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        out_ray = m.eval_child_ray(in_ray, 1, (0,0,0), (0,0,1))
        self.assertEquals(out_ray.direction, (1,0,-1))
        
        
class TestRay(unittest.TestCase):
    def test_termination(self):
        a = (1,2,3)
        b = (4,5,6)
        c = 2.3
        ray = ctracer.Ray(origin=a, direction=b)
        ray.length = c
        result = tuple(A+c*B for A,B in zip(a,b))
        shift = sum((a-b)**2 for a,b in zip(result, ray.termination))
        self.assertAlmostEquals(shift, 0.0)
        
        
class TestRayCollection(unittest.TestCase):
    def test_iteration(self):
        ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        rc = ctracer.RayCollection(10)
        for i in xrange(6):
            rc.add_ray(ray)
        self.assertEquals(rc.n_rays, 6)
        rays = [r for r in rc]
        self.assertEquals(len(rays), 6)
        
        
class AnOwner(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)
        
        
def sep(a,b):
    return sum((A-B)**2 for A,B in zip(a,b))
        
        
class TestTransform(unittest.TestCase):
    def test_trans(self):
        from enthought.tvtk.api import tvtk
        t = tvtk.Transform()
        t.rotate_x(10)
        t.rotate_y(15)
        t.rotate_z(20)
        t.translate(3,4,5)
        print t.matrix
        o = AnOwner(transform=t)
        fl = ctracer.FaceList(owner=o)
        fl.sync_transforms()
        print fl.transform.rotation
        print fl.transform.translation
        pt = (0.2, 0.4, 0.6)
        pt2 = t.transform_double_point(*pt)
        pt3 = ctracer.transform(fl.transform, pt)
        print pt, pt2, pt3
        self.assertEquals(sep(pt2,pt3), 0.0)
        
        
        
class TestTraceSegment(unittest.TestCase):
    def test_trace_segment(self):
        from raytrace import cfaces
        
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        
        rays = ctracer.RayCollection(10)
        in_ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        rays.add_ray(in_ray)
        
        self.assertEquals(rays.n_rays, 1)
        
        face_set = ctracer.FaceList()
        face_set.faces = [c]
        all_faces = [c]
        
        out_rays = ctracer.trace_segment(rays, [face_set], all_faces)
        
        self.assertEquals(out_rays.n_rays, 1)
        
        out_ray = out_rays[0]
        
        self.assertEquals(out_ray.direction, (1,0,-1))
        
if __name__=="__main__":
    unittest.main()
    

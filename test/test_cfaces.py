#!/usr/bin/env python

import pyximport
pyximport.install()

import sys
sys.path.append('..')
from raytrace import cfaces, ctracer
import unittest
from math import sqrt


class AnOwner(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class TestCircularFace(unittest.TestCase):        
    def test_params(self):
        c = cfaces.CircularFace()
        print c.params
        self.assertEquals(c.params, ['diameter', 'offset'])
        
    def test_update(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        c.update()
        self.assertEquals((c.diameter, c.offset), (o.diameter, o.offset))
        
    def test_intersection(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((0,0,-2), (0,0,2))
        self.assertEquals(dist, 2.0)
        
    def test_intersection2(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((-1,0,-1), (1,0,1))
        self.assertEquals(dist, sqrt(2.0))
        
    def test_face_set_intersection(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        c.idx = 7
        fl = ctracer.FaceList()
        fl.faces = [c]
        r = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1), length=10)
        idx = fl.intersect(r, 20)
        self.assertEquals(r.end_face_idx, 7)
        self.assertEquals(r.length, sqrt(2.0))
        
    def test_miss(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((6,0,-2),(6,0,2))
        import numpy
        self.assertEquals(dist, 0.0)
        
#    def test_eval_child(self):
#        o = AnOwner(diameter=5.5, offset=6.6)
#        c = cfaces.CircularFace(owner=o)
#        from raytrace import ctracer
#        in_ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
#        fl = ctracer.FaceList()
#        out_ray = c.eval_child_ray(in_ray, 1, (0,0,0), fl)
#        self.assertEquals(out_ray.direction, (1,0,-1))
    
if __name__=="__main__":
    unittest.main()
    
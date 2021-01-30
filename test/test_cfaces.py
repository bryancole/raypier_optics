#!/usr/bin/env python

#import pyximport
#pyximport.install()

import sys
# sys.path.append('..')
from raypier import cfaces, ctracer
import unittest
from math import sqrt


class AnOwner(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)

class TestCircularFace(unittest.TestCase):
    def test_params(self):
        c = cfaces.CircularFace()
        print(c.params)
        self.assertEqual(c.params, ['diameter', 'offset'])

    def test_update(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        c.update()
        self.assertEqual((c.diameter, c.offset), (o.diameter, o.offset))

    def test_intersection(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((0,0,-2), (0,0,2))
        self.assertEqual(dist, 2.0)

    def test_intersection2(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((-1,0,-1), (1,0,1))
        self.assertEqual(dist, sqrt(2.0))

    def test_face_set_intersection(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        c.idx = 7
        fl = ctracer.FaceList()
        fl.faces = [c]
        r = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1), length=10)
        idx = fl.intersect(r, 20)
        self.assertEqual(r.end_face_idx, 7)
        self.assertEqual(r.length, sqrt(2.0))

    def test_miss(self):
        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)
        dist = c.intersect((6,0,-2),(6,0,2))
        import numpy
        self.assertEqual(dist, 0.0)


class TestExtrudedFace(unittest.TestCase):
    def setUp(self):
        o = AnOwner()
        self.f = cfaces.ExtrudedPlanarFace(owner=o, z1=-1, z2=3,
                                            x1=-2.0,y1=2.0, x2=2.0, y2=-2.0)

    def test_intersection_1(self):
        dist = self.f.intersect((-5,0,0),(5,0,0))
        self.assertEqual(dist, 5.0)

    def test_miss_1(self):
        dist = self.f.intersect((-5,0,4),(5,0,4))
        self.assertEqual(dist, 0.0)

    def test_miss_2(self):
        dist = self.f.intersect((-5,0,-2),(5,0,-0.1))
        self.assertEqual(dist, 0.0)

    def test_miss_3(self):
        dist = self.f.intersect((-5,2.1,-1),(5,2.1,1))
        self.assertEqual(dist, 0.0)

if __name__=="__main__":
    unittest.main()

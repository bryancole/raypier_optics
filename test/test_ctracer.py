#!/usr/bin/env python

# import pyximport
# pyximport.install(language_level=3)

import sys
# sys.path.append('..')
from raypier import ctracer, cmaterials
import unittest
from math import sqrt
import random
import numpy


class TestRayStructure(unittest.TestCase):
    def test_sizeof_ray(self):
        self.assertEqual(ctracer.get_ray_size(), ctracer.ray_dtype.itemsize)


class TestVectorMathsFunctions(unittest.TestCase):
    def test_set_v(self):
        a = (1., 2., 3.)
        b = ctracer.py_set_v(a)
        self.assertEqual(a,b)

    def test_addvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.addvv(a,b)
        self.assertEqual(c, tuple(ai+bi for ai,bi in zip(a,b)))

    def test_invert(self):
        a = (1.2, 3.4, 5.6)
        b = ctracer.invert(a)
        self.assertEqual(a, tuple(-bi for bi in b))

    def test_subvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.subvv(a,b)
        self.assertEqual(c, tuple(ai-bi for ai,bi in zip(a,b)))

    def test_multvv(self):
        a = [1.,2.,3.]
        b = [2.,3.,4.]
        c = ctracer.multvv(a,b)
        self.assertEqual(c, tuple(ai*bi for ai,bi in zip(a,b)))

    def test_addvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.addvs(a,b)
        self.assertEqual(c, tuple(ai+b for ai in a))

    def test_subvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.subvs(a,b)
        self.assertEqual(c, tuple(ai-b for ai in a))

    def test_multvs(self):
        a = [1.,2.,3.]
        b = 4.56
        c = ctracer.multvs(a,b)
        self.assertEqual(c, tuple(ai*b for ai in a))

    def test_sep(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.sep(a,b)
        self.assertEqual(c, sqrt(sum((ai-bi)**2 for ai,bi in zip(a,b))))

    def test_dotprod(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.dotprod(a,b)
        self.assertEqual(c, sum(x*y for x,y in zip(a,b)))

    def test_cross(self):
        a = [1., 2., 3.]
        b = [4., 5., 6.]
        c = ctracer.cross(a,b)
        A = ctracer.dotprod(a,c)
        B = ctracer.dotprod(b,c)
        self.assertEqual(A, 0.0)
        self.assertEqual(B, 0.0)

    def test_map(self):
        a = [3.4, 5.6, 7.8]
        b = ctracer.mag(a)
        self.assertEqual(b, sqrt(sum(x**2 for x in a)))

    def test_norm(self):
        a = [1.1, 2.2, 3.3]
        b = ctracer.norm(a)
        print("B", b)
        self.assertAlmostEqual(1.0, sum(x**2 for x in b))
        self.assertAlmostEqual(0.0, ctracer.mag(ctracer.cross(a,b)))


class TestPECMaterial(unittest.TestCase):
    def test_eval_child(self):
        m = cmaterials.PECMaterial()
        in_ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        out_rays = ctracer.RayCollection(1)
        m.eval_child_ray(in_ray, 1, (0,0,0), (0,0,1), (0,1,0), out_rays)
        out_ray = out_rays[0]
        self.assertEqual(out_ray.direction, (1,0,-1))


class TestPolarisation(unittest.TestCase):
    def test_ellipticity_RHS(self):
        phase = (8.53 + 2.75j)
        in_ray = ctracer.Ray(origin=(-1,0,-1),
                             direction=(1,0,1),
                             E_vector=(1.,2.,3.),
                             E1_amp=phase*(1.0+0.0j)/numpy.sqrt(2),
                             E2_amp=phase*(0.0+1.0j)/numpy.sqrt(2))

        self.assertEqual(in_ray.ellipticity, 1.0)

    def test_ellipticity_LHS(self):
        phase = (1.23 + 3.45j)
        in_ray = ctracer.Ray(origin=(-1,0,-1),
                             direction=(1,0,1),
                             E_vector=(1.,2.,3.),
                             E1_amp=phase*(2.0+0.0j)/numpy.sqrt(2),
                             E2_amp=phase*(0.0-2.0j)/numpy.sqrt(2))

        self.assertEqual(in_ray.ellipticity, -1.0)

    def test_linear(self):
        in_ray = ctracer.Ray(origin=(-1,0,-1),
                             direction=(1,0,1),
                             E_vector=(1.,2.,3.),
                             E1_amp=(1.0+2.0j)/numpy.sqrt(2),
                             E2_amp=(2.0+4.0j)/numpy.sqrt(2))

        self.assertAlmostEqual(in_ray.ellipticity, 0.0)


class TestRay(unittest.TestCase):
    def test_termination(self):
        a = (1,2,3)
        b = (4,5,6)
        c = 2.3
        ray = ctracer.Ray(origin=a, direction=b)
        ray.length = c
        result = tuple(A+c*B for A,B in zip(a,b))
        shift = sum((a-b)**2 for a,b in zip(result, ray.termination))
        self.assertAlmostEqual(shift, 0.0)


class TestRayCollection(unittest.TestCase):
    def test_iteration(self):
        ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        rc = ctracer.RayCollection(10)
        for i in range(6):
            rc.add_ray(ray)
        self.assertEqual(rc.n_rays, 6)
        rays = [r for r in rc]
        self.assertEqual(len(rays), 6)

    def test_iteration2(self):
        ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        rc = ctracer.RayCollection(10)
        for i in range(6):
            rc.add_ray(ray)
        self.assertEqual(rc.n_rays, 6)
        itr = iter(rc)
        self.assertEqual(type(itr), ctracer.RayCollectionIterator)
        for i in range(6):
            r1 = next(itr)
            self.assertEqual(r1.origin, ray.origin)
            self.assertEqual(r1.direction, ray.direction)

    @staticmethod
    def make_ray():
        rnd = random.random
        r = ctracer.Ray(origin=(rnd(), rnd(), rnd()),
                        direction=(rnd(), rnd(), rnd()),
                        E1_amp = (rnd() + 1j*rnd()),
                        E2_amp = (rnd() + 1j*rnd()),
                        refractive_index = (rnd() + 1j*rnd())
                        )
        return r

    @staticmethod
    def compare(ray, row):
        a = tuple(ray.origin)==tuple(row['origin'])
        b = tuple(ray.direction)==tuple(row['direction'])
        return a and b

    def test_copy_to_array(self):
        rc = ctracer.RayCollection(10)
        rays = [self.make_ray() for i in range(7)]
        for ray in rays:
            rc.add_ray(ray)
        data = rc.copy_as_array()
        is_same = [self.compare(ray,row) for ray,row in zip(rays, data)]
        self.assertTrue(all(is_same))

    def test_from_array(self):
        a = numpy.empty(5, dtype=ctracer.ray_dtype)
        rc = ctracer.RayCollection.from_array(a)
        for i in range(len(a)):
            self.assertEqual(tuple(a[i]['origin']), tuple(rc[i].origin))
            self.assertEqual(tuple(a[i]['direction']), tuple(rc[i].direction))

    def test_properties(self):
        rc = ctracer.RayCollection(10)
        rays = [self.make_ray() for i in range(7)]
        for ray in rays:
            rc.add_ray(ray)
        data = rc.copy_as_array()
        self.assertTrue( numpy.alltrue( rc.E1_amp==data['E1_amp'] ) )
        self.assertTrue( numpy.alltrue( rc.E2_amp==data['E2_amp'] ) )
        self.assertTrue( numpy.alltrue( rc.refractive_index==data['refractive_index'] ) )


class AnOwner(object):
    def __init__(self, **kwds):
        self.__dict__.update(kwds)


def sep(a,b):
    return sum((A-B)**2 for A,B in zip(a,b))


class TestTransform(unittest.TestCase):
    """
    def test_trans(self):
        print "loading tvtk"
        from tvtk.api import tvtk
        print "complete"
        t = tvtk.Transform()
        t.rotate_x(10)
        t.rotate_y(15)
        t.rotate_z(20)
        t.translate(3,4,5)
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
    """


class TestTraceSegment(unittest.TestCase):
    def test_trace_segment(self):
        from raypier import cfaces

        o = AnOwner(diameter=5.5, offset=6.6)
        c = cfaces.CircularFace(owner=o)

        rays = ctracer.RayCollection(10)
        in_ray = ctracer.Ray(origin=(-1,0,-1), direction=(1,0,1))
        rays.add_ray(in_ray)

        self.assertEqual(rays.n_rays, 1)

        face_set = ctracer.FaceList()
        face_set.faces = [c]
        all_faces = [c]

        out_rays = ctracer.trace_segment(rays, [face_set], all_faces)

        self.assertEqual(out_rays.n_rays, 1)

        out_ray = out_rays[0]

        self.assertEqual(out_ray.direction, (1,0,-1))


class TestInterfaceMaterial(unittest.TestCase):
    def test_wavelengths(self):
        m = ctracer.InterfaceMaterial()
        a = [1.,2.,3.]
        m.wavelengths = numpy.array(a)
        b = list(m.wavelengths)

        self.assertEqual(a,b)

if __name__=="__main__":
    unittest.main()

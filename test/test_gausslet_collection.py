
import unittest
import numpy

from raypier.core.ctracer import GaussletCollection, gausslet_dtype, ray_dtype, RayCollection


class TestGaussletCollectionFromRays(unittest.TestCase):
    def test_from_rays(self):
        data = numpy.zeros(3, dtype=ray_dtype)
        data['origin'] = [[1.,2.,3.]]
        data['direction'] = [[4.,5.,6.]]
        data['length'] = [7.0,]
        
        gc = GaussletCollection.from_rays(data)
        rc = RayCollection.from_array(data)
        
        agc = gc[0]
        ar = rc[0]
        
        self.assertEqual(agc.base_ray.origin, ar.origin)
        self.assertEqual(agc.base_ray.direction, ar.direction)
        self.assertEqual(agc.base_ray.normal, ar.normal)
        self.assertEqual(agc.base_ray.length, ar.length)
        for i in range(6):
            paras = agc.parabasal_rays
            self.assertEqual(paras[i].origin, ar.origin)
            self.assertEqual(paras[i].direction, ar.direction)
            self.assertEqual(paras[i].length, ar.length)

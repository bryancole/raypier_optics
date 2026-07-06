
import unittest
import numpy as np

from tvtk.api import tvtk


class TestRotationEulerAngles(unittest.TestCase):
    def test_euler_conversions(self):
        o = 23.5
        e = 31.7
        r = 17.5
        t = tvtk.Transform()
        t.rotate_z(o)
        t.rotate_x(e)
        t.rotate_z(r)
        
        m = t.matrix.to_array()
        print(m)
        
        o2 = 180*np.arctan2(m[0,2], -m[1,2])/np.pi
        e2 = 180*np.arccos(m[2,2])/np.pi
        r2 = 180*np.arctan2(m[2,0],m[2,1])/np.pi
        
        print(o2,e2,r2)
        self.assertAlmostEqual(o, o2)
        self.assertAlmostEqual(e, e2)
        self.assertAlmostEqual(r, r2)
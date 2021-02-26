
import faulthandler
faulthandler.enable()

import unittest
import numpy

from matplotlib import pyplot as pp
from collections import Counter


from raypier.core.cdistortions import SimpleTestZernikeJ7, eval_nmk, zernike_R, zernike_Rprime, ZernikeDistortion
from raypier.core import cdistortions


class TestSimpleTestZernikeJ7(unittest.TestCase):
    def test_gradient_and_sag(self):
        dist = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=1.0)
        
        x_ = numpy.linspace(-10.,10.,200)
        y_ = numpy.linspace(-10.,10.,200)
        
        x,y = numpy.meshgrid(x_,y_)
        
        v = dist.z_offset_and_gradient(x.reshape(-1),y.reshape(-1))
        v = numpy.asarray(v)
        pp.imshow(v[:,2].reshape(x.shape))
        
        pp.show()


class TestEvalNMK(unittest.TestCase):
    def test_eval_nmk(self):
        n,m,k = eval_nmk(19)
        nmk = cdistortions.jnm_map
        
        ### Taken from  https://en.wikipedia.org/wiki/Zernike_polynomials
        ### ANSI standard J indexing
        ### table of (J,N,M)
        check = [(0,0,0),
                 (1,1,-1),
                 (2,1,1),
                 (3,2,-2),
                 (4,2,0),
                 (5,2,2),
                 (6,3,-3),
                 (7,3,-1),
                 (8,3,1),
                 (9,3,3),
                 (10,4,-4),
                 (11,4,-2),
                 (12,4,0),
                 (13,4,2),
                 (14,4,4),
                 (15,5,-5),
                 (16,5,-3),
                 (17,5,-1),
                 (18,5,1),
                 (19,5,3)
                  ]
        
        tab = [(i,n,m) for i,(n,m,k) in enumerate(nmk)]
        self.assertEqual(tab, check)
        
        tab = [(i,n,m) for i,(n,m,k) in enumerate(eval_nmk(i) for i in range(20))]
        self.assertEqual(tab, check)
        
        ###check there are no duplicate k-values for each (n, abs(m)) key.
        dups = {(n, abs(m)): k for n,m,k in nmk}
        c=Counter(dups.values())
        self.assertLessEqual(max(c.values()), 1)
        
        
class TestZernikeSeriesDistortion(unittest.TestCase):
    def test_eval_Z(self):
        Z = ZernikeDistortion(j7=1.0, unit_radius=10.0)
        Z2 = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=1.0)
        
        x_ = numpy.linspace(-5.,5.,4)
        x,y = numpy.meshgrid(x_,x_)
        x=x.ravel()
        y=y.ravel()
        a = Z.z_offset(x,y)
        b = Z2.z_offset(x,y)
        self.assertTrue(numpy.allclose(a, b))
        del Z
        
    def test_get_coefs(self):
        Z = ZernikeDistortion(j7=1.0, unit_radius=10.0)
        self.assertEqual(1.0, Z.get_coef(7))
        
    def test_eval_gradient(self):
        Z = ZernikeDistortion(j7=1.0, unit_radius=10.0)
        Z2 = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=1.0)
        
        x_ = numpy.linspace(-5.,5.,4)
        x,y = numpy.meshgrid(x_,x_)
        x=x.ravel()
        y=y.ravel()
        a = Z.z_offset_and_gradient(x,y)
        b = Z2.z_offset_and_gradient(x,y)
        self.assertTrue(numpy.allclose(a, b))
        del Z
        
    def test_gradient_at_origin(self):
        Z = ZernikeDistortion(j7=1.0, unit_radius=10.0)
        Z2 = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=1.0)
        x = y = numpy.array([0.0],'d')
        a = Z.z_offset_and_gradient(x,y)
        b = Z2.z_offset_and_gradient(x,y)
        print(a,b)
        self.assertTrue(numpy.allclose(a, b))
        
        
class TestEvalZernikeR(unittest.TestCase):
    def test_eval_zernike(self):
        r = 0.7071067811865476
        n,m,k = eval_nmk(7)
        k_max = max(eval_nmk(j)[2] for j in range(7))
        print("k_max=", k_max)
        workspace = numpy.zeros(k_max,'d')
        workspace[:] = float("nan")
        workspace[0] = 1.0
        print("k=",k, "n=",n, "m=",m)
        print("init=",workspace)
        R = zernike_R(r,k,n,m,workspace)
        print("R=",R)
        print("final=",workspace)
        
        R2 = 3*(r**3) - (2*r)
        self.assertAlmostEqual(R, R2)
        
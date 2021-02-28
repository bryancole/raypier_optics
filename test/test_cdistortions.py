
import faulthandler
faulthandler.enable()

import unittest
import numpy

from matplotlib import pyplot as pp
from collections import Counter


from raypier.core.cdistortions import SimpleTestZernikeJ7, eval_nmk, zernike_R, zernike_Rprime, \
        ZernikeDistortion, zernike_R_over_r
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
        
        
### Copied off Wikipedia
zernike_R_funcs = {
    (0,0): lambda r: 1.0,
    (1,1): lambda r: r,
    (2,0): lambda r: 2*(r**2) - 1,
    (2,2): lambda r: r**2,
    (3,1): lambda r: 3*(r**3) - (2*r),
    (3,3): lambda r: r**3,
    (4,0): lambda r: 6*(r**4) - 6*(r**2) + 1,
    (4,2): lambda r: 4*(r**4) - 3*(r**2),
    (4,4): lambda r: r**4,
    (5,1): lambda r: 10*(r**5) - 12*(r**3) + 3*r,
    (5,3): lambda r: 5*(r**5) - 4*(r**3),
    (5,5): lambda r: r**5,
    (6,0): lambda r: 20*(r**6) - 30*(r**4) + 12*(r**2) - 1,
    (6,2): lambda r: 15*(r**6) - 20*(r**4) + 6*(r**2),
    (6,4): lambda r: 6*(r**6) - 5*(r**4),
    (6,6): lambda r: r**6
    }


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
        
    def test_show_num_threads(self):
        print(cdistortions.num_threads)
        
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
        print(a)
        print(b)
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
        workspace = numpy.zeros((3,k_max),'d')
        workspace[0,:] = float("nan")
        workspace[0,0] = 1.0
        print("k=",k, "n=",n, "m=",m)
        print("init=",workspace)
        R = zernike_R(r,k,n,m,workspace)
        print("R=",R)
        print("final=",workspace)
        
        R2 = 3*(r**3) - (2*r)
        self.assertAlmostEqual(R, R2)
        
    def test_eval_zernike_R_all(self):
        r_all = numpy.linspace(0,1.0,20)
        k_max = 3*4 + 6 + 1
        workspace = numpy.zeros((3,k_max), 'd')
        for r in r_all: 
            workspace[:,:] = float("nan")
            for n in range(7):
                if 2*(n//2) == n:
                    start=0
                else:
                    start=1
                for m in range(start,n+1,2):
                    half_n = n//2
                    k=half_n*(half_n+1) + m
                    R = zernike_R(r,k,n,m,workspace)
                    R2 = zernike_R_funcs[(n,m)](r)
                    self.assertAlmostEqual(R,R2)
                
        
    def test_R_over_r(self):
        r=0.13
        n,m,k = eval_nmk(19)
        k_max = (6//2)*(6//2 + 1) + 6 + 1
        print(">> k_max=", k_max)
        k_max *= 2
        workspace = numpy.zeros((3,k_max),'d')
        workspace[:,:] = float("nan")
        #workspace[0,0] = 1.0

        for n in range(1,7):
            if 2*(n//2)==n:
                start=2
            else:
                start=1
            for m in range(start,n+1,2):
                half_n = int(n//2)
                k = half_n*(half_n+1) + abs(m)
                R = zernike_R(r,k,n,m,workspace)
                #Ra = zernike_R2(r,k,n,m,workspace)
                R2 = r*zernike_R_over_r(r,k,n,m,workspace)
                print("n,m,k=", (n,m,k), "R=",R, "R2=", R2)
                self.assertAlmostEqual(R,R2)
        
        
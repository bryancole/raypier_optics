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
        
if __name__=="__main__":
    unittest.main()
    

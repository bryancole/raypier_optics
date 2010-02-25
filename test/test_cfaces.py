#!/usr/bin/env python

import pyximport
pyximport.install()

import sys
sys.path.append('..')
from raytrace import cfaces
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
        
    
        
    
if __name__=="__main__":
    unittest.main()
    
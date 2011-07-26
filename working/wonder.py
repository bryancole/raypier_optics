#!/usr/bin/env python

#!/usr/bin/python
"""
A simple example
"""
import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import ConfocalRaySource
from raytrace.tracer import RayTraceModel
from raytrace.mirrors import RectWindow, RectMirror, PlanarWindow
from raytrace.lenses import PlanoConvexLens

import numpy


source = ConfocalRaySource(origin=(0,0,0),
 direction=(1,-2,0),
 working_dist = 80.,
 number=10,
 focus=(10,0,0),
 radius= 10,
 detail_resolution=5)

m1 = RectMirror(length=30,
 width = 30,
 thickness=5.0,
 centre=(0,0,0),
 direction=(1,0,1),
 n_inside = 1.333,
 name="m1")

'''                            
w1 = PlanarWindow(length=30,
 width = 15,
 radius=10,
 thickness=5.0,
 centre=(0,-20,0),
 direction=(0,1,1),
 n_inside = 1.333,
 name="w1")
'''
'''
l1 = PlanoConvexLens(diameter=25.4,
                thickness=6.0,
                centre=(0,-15,0),
                direction=(0,1,1),
                n_inside = 1.5,
                curvature = 20,
                name="l1")
'''
                
model = RayTraceModel(optics=[m1],
 sources=[source,])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

#import timeit
#t = timeit.Timer("model.update = True","from __main__ import model")
#ret = t.timeit(10)
#print "time:", ret



model.configure_traits()

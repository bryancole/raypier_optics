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
from raytrace.mirrors import PECMirror

from raytrace.cmaterials import PECMaterial

import numpy


source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=50,
                            rings=50,
                            detail_resolution=5,
                            theta=10.)
                            
m1 = PECMirror(diameter=25.4,
                thickness=5.0,
                centre=(0,-20,0),
                direction=(0,1,1),
                name="m1")

m2 = PECMirror(diameter=25.4,
                thickness=5.0,
                centre=(0,-20,-40),
                direction=(0,-1,-1),
                name="m2")
                
m3 = PECMirror(diameter=50,
                thickness=5.0,
                centre=(0,20,-40),
                direction=(0,-1,1),
                name="m3")
                
m4 = PECMirror(diameter=50,
                thickness=5.0,
                centre=(0,20,0),
                direction=(0,1,-1),
                name="m4")
                
model = RayTraceModel(optics=[m1,m2,m3,m4],
                    sources=[source,])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start            

model.update = True       

import timeit
t = timeit.Timer("model.update = True","from __main__ import model")
ret = t.timeit(10)
print("time:", ret)



model.configure_traits()

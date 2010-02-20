#!/usr/bin/env python

#!/usr/bin/python
"""
A simple example
"""
import sys
sys.path.append('..')
import pyximport
pyximport.install()

from raytrace.sources import ConfocalRaySource
from raytrace.tracer import RayTraceModel
from raytrace.mirrors import PECMirror

import numpy


source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=4,
                            rings=1,
                            detail_resolution=5,
                            theta=10.)
                            
m1 = PECMirror(diameter=25.4,
                thickness=5.0,
                centre=(0,-20,0),
                direction=(0,-1,0),
                name="m1")

print m1.faces.owner, m1
                
model = RayTraceModel(optics=[m1],
                    sources=[source,])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()
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

import numpy


source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=20,
                            detail_resolution=5,
                            theta=10.)

                
model = RayTraceModel(optics=[],
                    sources=[source,])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()
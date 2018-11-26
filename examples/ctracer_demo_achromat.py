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
from raytrace.achromats import EdmundOptic45805

import numpy


source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=4,
                            rings=1,
                            detail_resolution=5,
                            theta=10.)
                            
l1 = EdmundOptic45805(centre=(0,-20,0),
                      direction=(0,1,0),
                      )

                
model = RayTraceModel(optics=[l1,],
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
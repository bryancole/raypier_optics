#!/usr/bin/env python

#!/usr/bin/python
"""
A simple example
"""
import sys

import faulthandler
faulthandler.enable()

sys.path.append('..')
import pyximport
pyximport.install()

from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
from raytrace.achromats import EdmundOptic45805, Singlet
from raytrace.dispersion import NondispersiveCurve, NamedDispersionCurve


import numpy


# source = ConfocalRaySource(focus=(0,0,0),
#                             direction=(0,1,0),
#                             working_dist = 100.,
#                             number=4,
#                             rings=1,
#                             detail_resolution=5,
#                             theta=10.)

source = ParallelRaySource(direction=(0,1,0),
                           origin=(0,-50,0),
                           number=25,
                           radius=10,
                           rings=5)
                            
l1 = EdmundOptic45805(centre=(0,-20,0),
                      direction=(0,1,0),
                      )

# l1 = Singlet(center=(0,-20,0),
#              direction=(0,1,0),
#              CT=4.0,
#              curvature1 = 50,
#              curvature2 = -50,
#              dispersion = NondispersiveCurve(1.5),
#              dispersion_coating = NondispersiveCurve(1.25),
#              coating_thickness = 0.25
#              )

                
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

model.configure_traits(kind="live")

#!/usr/bin/python
"""
A simple example
"""
from raypier.sources import ParallelRaySource
from raypier.tracer import RayTraceModel
from raypier.troughs import TroughParabloid
import numpy

#
source = ParallelRaySource(origin=(0,20,0),
                            direction=(0,-1,1),
                            working_dist = 10.,
                            #number=20,
                            radius=1.)

#print source.InputDetailRays.origin.shape
                            
m1 = TroughParabloid(width = 20,
                       length = 100,
                       EFL = 5,
                       centre = (0,0,0))
                       
print(m1.intersect)
                
model = RayTraceModel(optics=[m1,],
                    sources=[source,],
                    probes=[],
                    recursion_limit=2)
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()

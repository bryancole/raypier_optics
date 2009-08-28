#!/usr/bin/python
"""
A simple example
"""
from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
#from raytrace.ellipsoids import Ellipsoid
#from raytrace.mirrors import PECMirror
#from raytrace.probes import PolarisationProbe
#from raytrace.lenses import PlanoConvexLens
from raytrace.involutes import CylindricalInvolute
import numpy

#p1 = PolarisationProbe(centre=(30,30,30),
#                        size=100.,
#                        direction=(1,1,1))
#

s1 = ParallelRaySource(origin=(100,100,30),
                            direction=(-.01,-1,0),
                            working_dist = 50.,
                            #number=20,
                            radius=14.)
"""                            
s2 = ParallelRaySource(origin=(x,y,z),
                            direction=(xd,yd,zd),
                            working_dist = 50.,
                            #number=20,
                            radius=0)
"""                            
#print source.InputDetailRays.origin.shape
                            
m1 = CylindricalInvolute(width = 60,
                       length = 500,
                       tube_radius = 5,
                       begin_angle = -90.,
                       end_angle = 500.)
                
                
model = RayTraceModel(optics=[m1,],
                    sources=[s1,],
                    probes=[])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()

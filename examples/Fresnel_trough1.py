#!/usr/bin/python
"""
A simple example of implementing a fresnel trough.  
"""
from raytrace.test_sources import PlaneRaySource
#from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
from raytrace.troughs import RectMirror
from raytrace.absorbers import RectAbsorber


import numpy


sun = PlaneRaySource(origin=(50,20,-100),
                            direction=(0,-1,0),
                            X_width = 100.,
                            X_resolution = 10.,
                            Y_width = 100.,
                            Y_resolution = 2.,
                            display = "wires")  
"""                            
sun = ParallelRaySource(origin=(50,20,-100),
                            direction=(0,-1,0),
                            radius = 20.,
                            display = "wires")                            
"""
mirror_width = 6.                          # 4 inches
centers = numpy.arange(0.1,100.,mirror_width)     # when all flat,they just touch.
tower_height = 50.                              # feet

# fake object to define array of objects later
blank = RectMirror(width = 20.,
                       length = 300.,
                       direction = (0,1,0),
                       centre = (0,0,0))

# find tilt for mirrors.
# (Functions applied to arrays of values)

sun_angle = numpy.pi/2                      #sun coming straight down


theta = numpy.arctan(tower_height/centers)
phi = theta/2 + sun_angle/2

print theta

# make an array of however many mirrors we have
Mirrors = list([blank] for i in centers)
    
for i,center in enumerate(centers):

     Mirrors[i] = RectMirror(width = mirror_width,
                       length = 300.,
                       direction = (center,numpy.tan(-phi[i])*center,0),
                       centre = (center,0,0))

Optics = Mirrors

sorber = RectAbsorber(width = mirror_width, length=300., direction = (0,-1,0),
                        centre = (0,tower_height,0))
Optics.append(sorber)
     
model = RayTraceModel(sources=[sun,], optics = Optics,)
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()

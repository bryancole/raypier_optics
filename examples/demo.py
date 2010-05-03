#!/usr/bin/python
"""
A simple example
"""
import sys
sys.path.append('..')
from raytrace.sources import ConfocalRaySource
from raytrace.tracer import RayTraceModel
from raytrace.ellipsoids import Ellipsoid
from raytrace.mirrors import PECMirror
#from raytrace.probes import PolarisationProbe
from raytrace.lenses import PlanoConvexLens
#from raytrace.absorbers import AbsorberDisk
import numpy

#p1 = PolarisationProbe(centre=(30,30,30),
#                        size=100.,
#                        direction=(1,1,1))

source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=20,
                            detail_resolution=5,
                            theta=10.)

#print source.InputDetailRays.origin.shape
                            
m1 = PECMirror(name="M1",
                centre=(0,-40,0),
                direction=(0,1,-1),
                diameter=30.,
                thickness=5.)
                
                
m4 = PECMirror(name="M4",
                centre=(0,40,0),
                direction=(0,-1,-1),
                diameter=30.,
                thickness=5.)
                
m2 = PECMirror(name="M2",
                centre=(0,-40,40),
                direction=(0,-1,1),
                diameter=15.,
                thickness=5.)
                
m2.orientation=-160
                
m3 = PECMirror(name="M3",
                centre=(0,40,40),
                direction=(0,1,1),
                diameter=15.,
                thickness=5.)
                
m3.orientation=20
                
m5 = PECMirror(name="M5",
                centre=(0,0,80),
                direction=(0,0,1),
                diameter=15.,
                thickness = 2.)
                
                
e1 = Ellipsoid(focus1=(0,-40,40),
                focus2=(0,0,80),
                size=80,
                X_bounds=(-30,0),
                Y_bounds=(-15,15),
                Z_bounds=(40-15.,40+15.))
                
e2 = Ellipsoid(focus1=(0,40,40),
                focus2=(0,0,80),
                size=80,
                X_bounds=(0,30),
                Y_bounds=(-15,15),
                Z_bounds=(40-15.,40+15.))      
                
l1 = PlanoConvexLens(centre=(0,-40,30),
                direction=(0,0,1))
                
model = RayTraceModel(optics=[m1,m2,m3,m4,m5,e1,e2, l1],
                    sources=[source,])
 
#model.trace_detail_async()
#import time
#start = time.clock()
#model.trace_all()
#end = time.clock()
#print "traced in", end - start                   

model.configure_traits()

#!/usr/bin/python
"""
A simple example
"""
from raytrace.sources import ConfocalRaySource
from raytrace.tracer import RayTraceModel
from raytrace.ellipsoids import Ellipsoid
from raytrace.mirrors import PECMirror

import numpy

source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            theta=10.)
                            
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
                
model = RayTraceModel(optics=[m1,m2,m3,m4,m5,e1,e2],
                    sources=[source,])
                    
model.configure_traits()

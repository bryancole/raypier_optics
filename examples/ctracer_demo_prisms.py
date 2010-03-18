#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
from raytrace.prisms import Rhomboid
from raytrace.cmaterials import PECMaterial, OpaqueMaterial

import numpy


source = ParallelRaySource(origin=(-7,-20,0),
                            direction=(0,1,0),
                            number=20,
                            radius=1.0,
                            scale_factor=0.1
                            )
                            
r1 = Rhomboid(name="rhomboid",
                #material=PECMaterial(),
                z_height_1=-5.0,
                z_height_2=5.0,
                height=7.0,
                width=14.0,
                slant=45.0,
                n_inside=1.5)
                            
                
model = RayTraceModel(optics=[r1,],
                    sources=[source,])
                    
model.configure_traits()
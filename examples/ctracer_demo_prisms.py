#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raypier.sources import ParallelRaySource
from raypier.tracer import RayTraceModel
from raypier.prisms import Rhomboid
from raypier.cmaterials import PECMaterial, OpaqueMaterial

import numpy


source = ParallelRaySource(origin=(-7,-20,0),
                            direction=(0,1,0),
                            number=20,
                            radius=1.0,
                            scale_factor=0.1
                            )
                            
r1 = Rhomboid(name="rhomboid",
                #material=OpaqueMaterial(),
                z_height_1=-5.0,
                z_height_2=5.0,
                height=7.0,
                width=14.0,
                slant=45.0,
                n_inside=1.5)
                            
                
model = RayTraceModel(optics=[r1,],
                    sources=[source,])
                    
model.configure_traits()
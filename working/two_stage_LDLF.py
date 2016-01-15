#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
from raytrace.prisms import LDLF
from raytrace.cmaterials import PECMaterial, OpaqueMaterial
from raytrace.constraints import BaseConstraint
from traits.api import Instance, on_trait_change, Tuple, Float
import numpy


source = ParallelRaySource(origin=(-2,30,0),
                            direction=(0,-1,0),
                            number=1,
                            radius=1.0,
                            scale_factor=0.1
                            )
			    
#linear dialectric filled light funnel                       
r2 = LDLF(#material=OpaqueMaterial(),
                centre = (0,5/numpy.sin(80.*numpy.pi/180.),0),
				z_height_1=-20.0,
                z_height_2=20.0,
                slat_width=10.0,
                ap_width=3+2*5*numpy.cos(80.*numpy.pi/180.),
                slant=75.0,
                n_inside=1.333) #1.333 = water

r1 = LDLF(#material=OpaqueMaterial(),
                centre = (0,0,0),
				z_height_1=-20.0,
                z_height_2=20.0,
                slat_width=5.0,	#
                ap_width=3.0,
                slant=80.0,
                n_inside=3) #1.5 = glass, 3 = something crazy

class Ap_match_Constraint(BaseConstraint):
   r1 = Instance(LDLF, ())
   r2 = Instance(LDLF, ())

   @on_trait_change("r1.update")
   def calc_r2_param(self):
		h = self.r1.slat_width/numpy.sin(self.r1.slant*numpy.pi/180.)
		self.r2.centre = (0,h,0)
		x = self.r1.ap_width + 2*self.r1.slat_width*numpy.cos(self.r1.slant*numpy.pi/180.)
		self.r2.ap_width = x

c = Ap_match_Constraint(r1=r1, r2=r2)

model = RayTraceModel(optics=[r1,r2], constraints=[c,],
                    sources=[source,])
                    
model.configure_traits()

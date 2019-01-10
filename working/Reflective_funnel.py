#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import ParallelRaySource, RectRaySource
from raytrace.tracer import RayTraceModel
from raytrace.prisms import Sheet
from raytrace.cmaterials import PECMaterial, OpaqueMaterial, TransparentMaterial
from raytrace.results import Ratio

import numpy

slat_width = 15
slant = 75
ex_ap = 3
ent_ap = ex_ap+2*slat_width*numpy.cos(slant*numpy.pi/180.)
height_offset = slat_width/numpy.sin(slant*numpy.pi/180.)


source = RectRaySource(origin=(0,30,15),
                            direction=(0,-1,0),
                            number=5,
                            length=10,
			    width = 3,
                            scale_factor=0.1
                            )

#sheet                 
r1 = Sheet(x1=ex_ap/2,
		y1=0,
		x2=ent_ap/2,
		y2=height_offset,
		z_height_1=0,
		z_height_2 = 30, 
		material = PECMaterial() )

r2 = Sheet(x1=-ex_ap/2,
		y1=0,
		x2=-ent_ap/2,
		y2=height_offset,
		z_height_1=0,
		z_height_2 = 30, 
		material = PECMaterial() )

en_ap = Sheet(x1=ent_ap/2,
		y1=height_offset,
		x2=-ent_ap/2,
		y2=height_offset,
		z_height_1=0,
		z_height_2 = 30, 
		material = TransparentMaterial() )

exit_ap = Sheet(x1=ex_ap/2,
		y1=0,
		x2=-ex_ap/2,
		y2=0,
		z_height_1=0,
		z_height_2 = 30, 
		material = TransparentMaterial() )

                           
ratio = Ratio(denominator=en_ap.faces[0],
              nominator=exit_ap.faces[0])

model = RayTraceModel(optics=[r1,r2,en_ap,exit_ap],
                    sources=[source,], results=[ratio])
                    
model.configure_traits()

print(model.all_faces)
ap_ratio = ent_ap/ex_ap
print("ap ratio: ", ap_ratio)
print("effective concentration: ", ap_ratio*model.results[0].result)	#still not right, circular source biases ratio.

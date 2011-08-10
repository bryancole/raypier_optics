#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import RectRaySource
from raytrace.tracer import RayTraceModel
from raytrace.dielectrictroughs import LDLF
from raytrace.cmaterials import PECMaterial, OpaqueMaterial
from raytrace.results import Ratio

import numpy

slat_width = 20
slant = 79
ex_ap = 3
ent_ap = ex_ap+2*slat_width*numpy.cos(slant*numpy.pi/180.)
height_offset = slat_width/numpy.sin(slant*numpy.pi/180.)


source = RectRaySource(    origin = (0.,height_offset+.1,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap,
                rings = 3,
                theta = -10,
                randomness = True)

#linear dialectric filled light funnel                       
r1 = LDLF(#material=OpaqueMaterial(),
                z_height_1=-20.0,
                z_height_2=20.0,
                slat_width=slat_width,
                ap_width=ex_ap,
                slant=slant,
                n_inside=1.333) #1.333 = water


                           
ratio = Ratio(denominator=r1.faces[2],
              nominator=r1.faces[0])

model = RayTraceModel(optics=[r1,],
                    sources=[source,], results=[ratio,],recursion_limit=5)
                    
#model.configure_traits()
print model.all_faces

ap_ratio = ent_ap/ex_ap
print "ap ratio: ", ap_ratio
print "effective concentration: ", ap_ratio*model.results[0].result	#still not right, circular source biases ratio.

model.configure_traits()

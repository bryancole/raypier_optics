#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import RectRaySource
from raytrace.tracer import RayTraceModel
from raytrace.dielectrictroughs import LDLF
from raytrace.cmaterials import PECMaterial, OpaqueMaterial
from raytrace.results import Ratio, RayPaths
from raytrace.render_utils import trace_pretty
import numpy

slat_width = 20
slant = 79
ex_ap = 3
ent_ap = ex_ap+2*slat_width*numpy.cos(slant*numpy.pi/180.)
height_offset = slat_width/numpy.sin(slant*numpy.pi/180.)


s1 = RectRaySource(    origin = (0.,2*height_offset,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap+10,
                theta = -15,
                randomness = True)

s2 = RectRaySource(    origin = (0.,2*height_offset,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap+10,
                theta = -5,
                randomness = True)
s3 = RectRaySource(    origin = (0.,2*height_offset,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap+10,
                theta = 0,
                randomness = True)
s4 = RectRaySource(    origin = (0.,2*height_offset,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap+10,
                theta = 5,
                randomness = True)
s5 = RectRaySource(    origin = (0.,2*height_offset,0.),
                direction = (0,-1,0),
                number = 20,		#total rays is n^2
                length = ent_ap,
                width = ent_ap+10,
                theta = 15,
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
raypaths = RayPaths()

model = RayTraceModel(optics=[r1,],
                    sources=[s1,s2,s3,s4,s5], results=[ratio,raypaths],recursion_limit=5)

entry_ap = r1.faces.faces[2]
newsource = trace_pretty(model,1,[entry_ap])                    

model.sources = [newsource,]
#model.configure_traits()
#print model.all_faces

#ap_ratio = ent_ap/ex_ap
#print "ap ratio: ", ap_ratio
#print "effective concentration: ", ap_ratio*model.results[0].result	#still not right, circular source biases ratio.

model.configure_traits()

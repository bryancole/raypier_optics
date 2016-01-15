#!/usr/bin/env python

import sys
sys.path.append('..')
#import pyximport
#pyximport.install()

from raytrace.sources import ParallelRaySource
from raytrace.tracer import RayTraceModel
from raytrace.prisms import Extrusion
from raytrace.cmaterials import PECMaterial, OpaqueMaterial
from raytrace.results import Ratio

#mass copy from other file, this is not all needed:
from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Bool
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, izip, islice, tee
from raytrace.bases import Optic, Traceable
from raytrace.cfaces import PolygonFace, ExtrudedPlanarFace
from raytrace.ctracer import FaceList
##



source = ParallelRaySource(origin=(0,30,0),
                            direction=(0,-1,0),
                            number=2000,
                            radius=2,
			    rings = 1,
                            scale_factor=0.1)

class approx_parabola(Extrusion):
    name = "A polygon approximation of a parabolic mirror"
    focal_length = Float #Focal distance of parabola
    concave = Bool(True) #concave or convex surface exposed
    left_most = Float(-10)
    right_most = Float(10)
    resolution = Float(.01)
    Extrusion.trace_ends = Bool(False)

    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('trace_ends'),
                       #Item('n_inside'),
                       Item('focal_length'),
                       Item('length')))

    def parabola(self, x):
	pt = 4*x*x/self.focal_length
    	return pt                       

    @on_trait_change("focal_length")
    def config_profile(self):
        fl = self.focal_length
	neg = self.left_most
	pos = self.right_most
	res = self.resolution

	#generate array of points by sampling the parabolic equation
	profile = []
	abscissa = numpy.arange(neg,pos,res)
	for i,x in enumerate(abscissa):
		profile.append((x,self.parabola(x)))
	if self.concave:
		xtr_pts = [(pos,0),(0,-1),(neg,0)]
		profile[len(profile):] = xtr_pts

	self.profile=profile
	#print profile

#polygon approximation of parabolic trough mirror                 
para = approx_parabola(material=PECMaterial(),
                z_height_1=-20.0,
                z_height_2=20.0,
                focal_length=10)


          
print 0
model0 = RayTraceModel(optics=[para,], sources=[source,])
print 1
model1 = RayTraceModel(optics=[para,], sources=[source,])
print 2
model2 = RayTraceModel(optics=[para,], sources=[source,])
print 3
                    
model1.configure_traits()
print model.all_faces.size


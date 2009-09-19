"""
A demonstration of how to use a Results subclass
"""


from raytrace.sources import ConfocalRaySource
from raytrace.tracer import RayTraceModel
from raytrace.mirrors import PECMirror
from raytrace.lenses import PlanoConvexLens
from raytrace.results import Total_Efficency, Ratio
from raytrace.apertures import RectAperture

import numpy

source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,-1,0),
                            working_dist = 100.,
                            number=5,
                            detail_resolution=5,
                            theta=10.)

#print source.InputDetailRays.origin.shape
                            
m1 = PECMirror(name="M1",
                centre=(0,-40,0),
                direction=(0,1,-1),
                diameter=30.,
                thickness=5.)
                
aligned = RectAperture(name="entry",
                    centre = (0,5,15),
                    direction = (0,1,0),
                    length = 30.,
                    width = 30.)
                    
antialigned = RectAperture(name="reversed",
                    centre = (0,-5,15),
                    direction = (0,-1,0),
                    length = 30.,
                    width = 30.)
                
l1 = PlanoConvexLens(centre=(0,-40,30),
                direction=(0,0,1))
                
ratio = Total_Efficency(Target=m1.faces[0],
              Aperture=aligned.faces[0])

model = RayTraceModel(optics=[m1,aligned,antialigned],
                    sources=[source,],
                    results=[ratio])               

model.configure_traits()
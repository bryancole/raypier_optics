"""
A demonstration of how to use a Results subclass
"""
import sys
sys.path.append('..')

from raypier.sources import ConfocalRaySource
from raypier.tracer import RayTraceModel
from raypier.mirrors import PECMirror
from raypier.lenses import PlanoConvexLens
from raypier.results import Ratio

import numpy

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
    
                
l1 = PlanoConvexLens(centre=(0,-40,30),
                direction=(0,0,1))

ratio = Ratio(denominator=m1.faces[0],
              nominator=l1.faces[1])
                
model = RayTraceModel(optics=[m1,l1],
                    sources=[source,],
                    results=[ratio])               

model.configure_traits()
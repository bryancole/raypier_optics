#!/usr/bin/env python

#!/usr/bin/python
"""
A simple example
"""
import sys
sys.path.append('..')

from raypier.sources import ConfocalRaySource
from raypier.tracer import RayTraceModel
from raypier.lenses import PlanoConvexLens


source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=50,
                            rings=16,
                            detail_resolution=5,
                            theta=10.)
                            
l1 = PlanoConvexLens(diameter=25.4,
                thickness=6.0,
                centre=(0,-20,0),
                direction=(0,1,0),
                n_inside = 1.5,
                curvature = 20,
                name="l1")

                
model = RayTraceModel(optics=[l1,],
                    sources=[source,])
 

model.configure_traits()
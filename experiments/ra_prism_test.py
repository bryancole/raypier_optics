

import sys
import numpy
#from matplotlib import pyplot as pp
print(numpy.__version__)

from traits.api import Bool, Float, observe, Range
from traitsui.api import VGroup, View, Item

from raypier.tracer import RayTraceModel, Traceable

from raypier.prisms import RightAnglePrism

from raypier.sources import ParallelRaySource, SingleRaySource

beam_offset=0.0
source = ParallelRaySource(origin=(-50.0,-beam_offset,0.0),
                           direction=(1.,0.,0.), #(0.,-0.008,-1.),
                           radius=1.0,
                           scale_factor=0.5,
                           number=8,
                           rings=5,
                           E_vector=(0.,0.,1.),
                           E1_amp = 1.0,
                           E2_amp = 0.0,
                           display='wires'
                           )

side = 10.0
hypotenuse = numpy.sqrt(2*(side**2))
M0 = RightAnglePrism(centre=(0.0, 0.0, 0.0),
                     #direction=(0,0,-1),
                     direction=(-1,0,0),
                     rotation=-90.0,
                     z_height_1=-side/2,
                     z_height_2=side/2,
                     width=hypotenuse,
                     height=hypotenuse/2.,
                     n_inside=1.51 #N-BK7
                     )

model = RayTraceModel(sources=[source],
                      optics=[M0])
model.recursion_limit = 10
model.configure_traits()

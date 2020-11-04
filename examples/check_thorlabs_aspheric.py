
from raytrace.aspherics import Thorlabs355440
from raytrace.tracer import RayTraceModel
from raytrace.sources import ConfocalRaySource

from raytrace.find_focus import find_ray_focus

import numpy

lens = Thorlabs355440(centre=(0,0,0),
                    direction=(1,0,0))

source = ConfocalRaySource(focus = (-7.09,0,0),
                           direction = (1,0,0),
                           working_dist= 10.0,
                           theta=12.0)
 
model = RayTraceModel(optics=[lens],
                      sources=[source])

x1 = numpy.linspace(2.0,8.0,500)
x2 = []

for x in x1:
    source.focus = (-x,0,0)
    source.theta = 6.0 * numpy.arctan2(2.,x)/0.2
    model.trace_all()
    last_rays = source.TracedRays[-1]
    focus = find_ray_focus(last_rays)
    x2.append( focus[0] - lens.CT )
    print(x, source.theta,  x2[-1] , focus)

from matplotlib import pyplot as pp
pp.plot(x1,x2,'-')

pp.show()

model.configure_traits()


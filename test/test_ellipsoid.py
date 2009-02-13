from raytrace.tracer import RayTraceModel
from raytrace.ellipsoids import Ellipsoid
from raytrace.sources import ParallelRaySource, ConfocalRaySource

ellipse = Ellipsoid(focus1=(-50,0,0),
                    focus2=(0,50,0),
                    Z_bounds=(-25,25),
                    size=100.)

#source = ParallelRaySource(origin=(-60,0,12),
#                           direction=(1,0,0),
#                           radius=10.)
source = ConfocalRaySource(focus=(-50,0,0),
                           direction=(1,0,0),
                           theta=15.,
                           working_dist=10.)

model = RayTraceModel(optics=[ellipse,],
                      sources=[source,])

model.configure_traits()

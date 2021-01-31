from raypier.tracer import RayTraceModel
from raypier.ellipsoids import Ellipsoid
from raypier.sources import ParallelRaySource, ConfocalRaySource

e = Ellipsoid(focus1=(-50,0,0),
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

model = RayTraceModel(optics=[e,],
                      sources=[source,])

#model.update = True
#
#print source.traced_rays[0]

model.configure_traits()

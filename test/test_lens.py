from raytrace.tracer import RayTraceModel
from raytrace.lenses import PlanoConvexLens
from raytrace.sources import ConfocalRaySource

lens = PlanoConvexLens(centre=(20,0,0),
                       direction=(1,0,0),
                       curvature=25.0,
                       diameter=25.4,
                       CT=6.0)

source = ConfocalRaySource(focus = (0,0,0),
                           direction = (1,0,0),
                           working_dist= 100.0,
                           theta=15.0)

model = RayTraceModel(optics=[lens],
                      sources=[source])

model.configure_traits()

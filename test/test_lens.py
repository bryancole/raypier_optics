from raypier.tracer import RayTraceModel
from raypier.lenses import PlanoConvexLens
from raypier.sources import ConfocalRaySource

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

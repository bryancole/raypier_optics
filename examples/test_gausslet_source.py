

from raytrace.gausslet_sources import SingleGaussletSource
from raytrace.lenses import PlanoConvexLens

from raytrace.tracer import RayTraceModel

lens = PlanoConvexLens(centre=(0,0,50),
                       n_inside=1.5,
                       curvature=15.0)
src = SingleGaussletSource(beam_waist=1.0, max_ray_len=55.0)

model = RayTraceModel(sources=[src], optics=[lens])
model.configure_traits()
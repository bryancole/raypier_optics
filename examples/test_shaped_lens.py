

from raytrace.tracer import RayTraceModel
from raytrace.lenses import ShapedPlanoSphericLens
from raytrace.sources import ConfocalRayFieldSource
from raytrace.mirrors import SphericalMirrorWithHole


lens = ShapedPlanoSphericLens(centre=(0,0,150), direction=(0,0,1),
                              diamter=50.0,
                              curvature=25, CT=6, n_inside=1.5)

m = SphericalMirrorWithHole(centre=(0,0,200))

src = ConfocalRayFieldSource()


model = RayTraceModel(optics=[lens,m], sources=[src])

model.configure_traits()

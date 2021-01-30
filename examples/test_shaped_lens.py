

from raypier.tracer import RayTraceModel
from raypier.lenses import ShapedPlanoSphericLens
from raypier.sources import ConfocalRayFieldSource
from raypier.mirrors import SphericalMirrorWithHole


lens = ShapedPlanoSphericLens(centre=(0,0,150), direction=(0,0,1),
                              diamter=50.0,
                              curvature=25, CT=6, n_inside=1.5)

m = SphericalMirrorWithHole(centre=(0,0,200))

src = ConfocalRayFieldSource()


model = RayTraceModel(optics=[lens,m], sources=[src])

model.configure_traits()

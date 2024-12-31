
from raypier.tracer import RayTraceModel
from raypier.bezier2d import BSplineSinglePatch

from raypier.sources import ParallelRaySource
from raypier.gausslet_sources import CollimatedGaussletSource


patch = BSplineSinglePatch()


if True:
    src = ParallelRaySource(origin=(0,0,-50.0),
                            direction=(0,0,1),
                            radius = 3.0,
                            rings=2,
                            display="pipes",
                            show_normals=True,
                            opacity=0.4)
else:
    src = CollimatedGaussletSource(origin=(0,0,-50.0),
                            direction=(0,0,1),
                            radius = 3.0,
                            resolution=0.5,
                            display="pipes",
                            show_normals=True,
                            opacity=0.4)

model = RayTraceModel(optics=[patch], sources=[src], recursion_limit=5)
model.configure_traits()
 
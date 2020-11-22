
from raytrace.tracer import RayTraceModel
from raytrace.sources import HexagonalRayFieldSource
from raytrace.lenses import PlanoConvexLens
from raytrace.fields import EFieldPlane

lens = PlanoConvexLens(centre=(0,0,20),
                       direction=(0,0,1),
                       diameter=25.0,
                       thickness=6.0,
                       curvature=40.0,
                       n_inside=1.5)

src = HexagonalRayFieldSource(spacing=2.0, direction=(0,0,1))

probe = EFieldPlane(source=src,
                    centre=(0,0,70),
                    direction=(0,0,1),
                    width=30,
                    size=30)

model = RayTraceModel(sources=[src], optics=[lens],
                      probes=[probe])

model.configure_traits()

for item in src.iter_neighbours():
    print(item)
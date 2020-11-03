import pyximport
pyximport.install(language_level=3)

from raytrace.tracer import RayTraceModel
from raytrace.lenses import PlanoConicLens, PlanoConvexLens, BiConicLens, AsphericLens
from raytrace.sources import ConfocalRaySource

lens = PlanoConicLens(centre=(20,0,0),
                       direction=(1,0,0),
                       curvature=25.0,
                       diameter=25.4,
                       conic_const = 0.0,
                       CT=6.0)
# 
# lens = AsphericLens(centre=(20,0,0),
#                     direction=(1,0,0),
#                     CT=12.0,
#                     )

lens = BiConicLens(centre=(20,0,0),
                       direction=(1,0,0),
                       curvature=25.0,
                       curvature2=-25.0,
                       diameter=25.4,
                       conic_const = 0.0,
                       conic_const2 = 0.0,
                       CT=12.0)

# lens2 = PlanoConvexLens(centre=(20,0,0),
#                        direction=(1,0,0),
#                        curvature=25.0,
#                        diameter=25.4,
#                        conic_const = 0.0,
#                        CT=6.0)

#print("params", lens.params)
for f in lens.faces.faces:
    print(f.params)
 
source = ConfocalRaySource(focus = (0,0,0),
                           direction = (1,0,0),
                           working_dist= 100.0,
                           theta=15.0)
 
model = RayTraceModel(optics=[lens],
                      sources=[source])

model.configure_traits()

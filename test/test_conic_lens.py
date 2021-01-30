import pyximport
pyximport.install(language_level=3)

from raypier.tracer import RayTraceModel
from raypier.lenses import PlanoConicLens, PlanoConvexLens, BiConicLens, AsphericLens
from raypier.sources import ConfocalRaySource

lens = PlanoConicLens(centre=(20,0,0),
                       direction=(1,0,0),
                       curvature=25.0,
                       diameter=25.4,
                       conic_const = 0.0,
                       CT=6.0)
# 
# lens = AsphericLens(centre=(20,0,0),
#                     direction=(1,0,0),
#                     CT=3.8*5,
#                     diameter=4.7*5,
#                     A_curvature = -5*3.200655,
#                     A_conic = -4.321649,
#                     A4 = -5.521153E-3/(5**3),
#                     A6 = 1.981378E-3/(5**5),
#                     A8 = -4.782553E-4/(5**7),
#                     A10 = 7.328134E-5/(5**9),
#                     B_curvature = 5*3.200655,
#                     B_conic = -4.321649,
#                     B4 = 5.521153E-3/(5**3),
#                     B6 = -1.981378E-3/(5**5),
#                     B8 = 4.782553E-4/(5**7),
#                     B10 = -7.328134E-5/(5**9),
#                     )
# 
# lens = BiConicLens(centre=(20,0,0),
#                        direction=(1,0,0),
#                        curvature=25.0,
#                        curvature2=-25.0,
#                        diameter=25.4,
#                        conic_const = 0.0,
#                        conic_const2 = 0.0,
#                        CT=12.0)

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

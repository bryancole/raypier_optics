from raypier.tracer import RayTraceModel
from raypier.parabolics import OffAxisParabloid
from raypier.sources import ConfocalRaySource

h = 25.4

oap = OffAxisParabloid(diameter=50.4,
                       height = h,
                       EFL = 50.4)

source = ConfocalRaySource(focus = (0,0,0),
                           direction = (1,0,0),
                           working_dist= 10.0,
                           number=8,
                           rings=2,
                           theta=15.0,
                           export_pipes=True)

model = RayTraceModel(optics=[oap],
                      sources=[source])

model.configure_traits()
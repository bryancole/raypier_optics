from raytrace.tracer import RayTraceModel
from raytrace.parabolics import OffAxisParabloid
from raytrace.sources import ConfocalRaySource

h = 25.4

oap = OffAxisParabloid(diameter=50.4,
                       height = h,
                       EFL = 50.4)

source = ConfocalRaySource(focus = (0,0,h),
                           direction = (1,0,0),
                           working_dist= 10.0,
                           theta=15.0)

model = RayTraceModel(optics=[oap],
                      sources=[source])

model.configure_traits()
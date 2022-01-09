
from raypier.tracer import RayTraceModel
from raypier.bezier2d import BSplinePatch


patch = BSplinePatch()

model = RayTraceModel(optics=[patch])
model.configure_traits()
 
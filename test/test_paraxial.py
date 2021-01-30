from raypier.tracer import RayTraceModel
from raypier.paraxial import ParaxialLens
from raypier.sources import ConfocalRaySource

lens = ParaxialLens(diameter = 50.0,
                    focal_length = 50.0,
                    centre = (0,0,0),
                    direction = (1,0,0))

lens2 = ParaxialLens(diameter = 50.0,
                    focal_length = 50.0,
                    centre = (50.0,0,0),
                    direction = (1,0,0))

source = ConfocalRaySource(focus = (-50.0,0,0),
                           direction = (1,0,0),
                           working_dist= 10.0,
                           theta=15.0)

model = RayTraceModel(optics=[lens,lens2],
                      sources=[source])

model.configure_traits()

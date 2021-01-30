

from raypier.apertures import CircularAperture

from raypier.tracer import RayTraceModel

from raypier.sources import ConfocalRayFieldSource


src = ConfocalRayFieldSource(centre=(0,0,0),
                             direction=(0,0,1),
                             working_dist=50,
                             angle=5.0,
                             angle_step=1.0)

ap = CircularAperture(centre=(0,0,10),
                      direction=(0,0,1),
                      outer_diameter=25.0,
                      inner_diameter=15.0,
                      hole_diameter=5,
                      edge_width=2.0)

model = RayTraceModel(optics=[ap],
                      sources=[src])

model.configure_traits()

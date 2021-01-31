from raypier.lenses import PlanoConvexLens
from raypier.sources import GaussianBeamRaySource
from raypier.tracer import RayTraceModel

m1 = PlanoConvexLens(name="L1",
                centre=(0,-45, 0),
                direction=(0.,-1,0), #(0,y_dir,-small_height),
                diameter=10.,
                CT=4.,
                curvature=15.0,
                n_inside=1.6)

source = GaussianBeamRaySource(origin=(0,0,0),
                           direction=(0,-1,0),
                           beam_waist=1000.0, #in microns
                           number=20,
                           wavelength=100)

model = RayTraceModel(sources=[source],
                      optics=[m1])

model.configure_traits()
#model.update = True
#print source.traced_rays
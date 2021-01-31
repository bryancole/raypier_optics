from raypier.mirrors import PECMirror
from raypier.sources import ParallelRaySource
from raypier.tracer import RayTraceModel

m1 = PECMirror(name="M1",
                centre=(0,-45, 0),
                direction=(0.,-1,1), #(0,y_dir,-small_height),
                diameter=40.,
                thickness=4.)

source = ParallelRaySource(origin=(0,0,0),
                           direction=(0,-1,0),
                           radius=12,
                           number=20)

model = RayTraceModel(sources=[source],
                      optics=[m1])

model.configure_traits()
#model.update = True
#print source.traced_rays
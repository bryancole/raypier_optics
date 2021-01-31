from raypier.mirrors import PECMirror
from raypier.sources import ConfocalRaySource
from raypier.tracer import RayTraceModel

m1 = PECMirror(name="M1",
                centre=(0,-45, 0),
                direction=(0.,1,-1), #(0,y_dir,-small_height),
                diameter=40.,
                thickness=4.)

source = ConfocalRaySource(focus=(0,0,0),
                            direction=(0,1,0),
                            working_dist = 100.,
                            number=20,
                            detail_resolution=5,
                            theta=15.00, #fringe ray is at 15.22 deg
                            scale_factor=0.2)

model = RayTraceModel(sources=[source],
                      optics=[m1])

model.configure_traits()
#model.update = True
#print source.traced_rays


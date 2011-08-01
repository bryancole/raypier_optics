from raytrace.splines import Extruded_interpolant
from raytrace.sources import RectRaySource
from raytrace.tracer import RayTraceModel
from raytrace.cmaterials import OpaqueMaterial, PECMaterial
import numpy as np

x = np.linspace(1, 30, 100)
# dont use form like this:
#x = np.array([1,2,3,4,5,6,7,8,9,10])
# because these are used as integers and decimals are not used!
y = -30/x 

#give profile as list of x values then y values.  this profile is evaluated by splprep
#the smaller smoothness is, the more accurate the spline is, but the more faces it takes.
#try to make s as small as possible without making too many faces. (20 faces isn't unreasonable)
test = Extruded_interpolant(profile = [x,y], smoothness = .5, z_height_1 = -30, z_height_2=30, \
		material = PECMaterial(), trace_ends= False, trace_top = False )

source = RectRaySource(origin=(5,50,0),
                            direction=(0,-.5,0),
                            working_dist = 100.,
                            number=1,
                            length = 10,
            			    width = 10,
            			    randomness = False,
                            scale_factor=0.2)

model = RayTraceModel(sources=[source],
                      optics=[test], recursion_limit=2)
model.configure_traits()


from raypier.splines import Extruded_bezier, b_spline_to_bezier_series
from raypier.sources import RectRaySource
from raypier.tracer import RayTraceModel
from raypier.core.cmaterials import OpaqueMaterial, PECMaterial,DielectricMaterial
from raypier.results import RayPaths
import numpy as np

from scipy.interpolate import splprep, splev

x = np.linspace(1, 30, 100)
# dont use form like this:
#x = np.array([1,2,3,4,5,6,7,8,9,10])
# because these are used as integers and decimals are not used!
y = 30/x 


tck, uout = splprep([x,y], s=.0005, k=3, per=False)
ctrl_pts = b_spline_to_bezier_series(tck)


#give profile as list of x values then y values.  this profile is evaluated by splprep
#the smaller smoothness is, the more accurate the spline is, but the more faces it takes.
#try to make s as small as possible without making too many faces. (20 faces isn't unreasonable)
test = Extruded_bezier(control_points = ctrl_pts, z_height_1 = -30, z_height_2=29, \
		material = DielectricMaterial(), trace_ends= True, trace_top = True, invert_normal=True)
'''
flattest = Extruded_bezier(control_points = np.array([[[0,0],[10,0],[20,0],[30,0]]]), z_height_1 = -30, z_height_2=30, \
		material = DielectricMaterial(), trace_ends= True, trace_top = False, invert_normal=True )
'''

source = RectRaySource(origin=(10,30,0),
                            direction=(.1,-1,0),
                            working_dist = 100.,
                            number=1,
                            length = 0,
            			    width = 0,
            			    randomness = False,
                            scale_factor=0.2)
'''

source = RectRaySource(origin=(10.,10.,0),
                            direction=(.2,-1,0),
                            working_dist = 100.,
                            number=3,
                            length = 12,
            			    width = 12,
            			    randomness = False,
                            scale_factor=0.2)
'''
raypaths = RayPaths()

model = RayTraceModel(sources=[source],
                      optics=[test], results = [raypaths], recursion_limit=10)

#test.invert_normal = True
hits = model.sources[0].traced_rays
res = model.results[0].result
model.configure_traits()



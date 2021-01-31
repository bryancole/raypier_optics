from raypier.splines import Extruded_interpolant
from raypier.sources import RectRaySource
from raypier.tracer import RayTraceModel
from raypier.core.cmaterials import OpaqueMaterial, PECMaterial, DielectricMaterial
from raypier.core.cfaces import ExtrudedPlanarFace
from raypier.dielectrictroughs import LDLF

import numpy as np

t = np.linspace(0, 2*np.pi, 100)
# dont use form like this:
#x = np.array([1,2,3,4,5,6,7,8,9,10])
# because these are used as integers and decimals are not used!
# range() also creates a list of integers, so use numpy.linspace
#y = -x*x+30*30*x

xs = 30*np.sin(t)
ys = 30*np.cos(t)

for s,x in enumerate(xs):
    if x<0:
        xs[s]=x/2. 
xs = xs+50
#give profile as list of x values then y values.  this profile is evaluated by splprep
#the smaller smoothness is, the more accurate the spline is, but the more faces it takes.
#try to make s as small as possible without making too many faces. (20 faces isn't unreasonable)

test = Extruded_interpolant(profile = [xs,ys], smoothness = 0, z_height_1 = -30, z_height_2=30, \
		material = DielectricMaterial(), n_inside=1.5, trace_ends= True,trace_top = False )

source = RectRaySource(origin=(50,0,0),
                            direction=(-1,0,0),
                            working_dist = 100.,
                            number=10,
                            length = 10,
			    width = 10,
			    randomness = True,
                            theta = 5.,
                            scale_factor=0.2)

'''
r1 = LDLF(material=PECMaterial(),
                z_height_1=-20.0,
                z_height_2=20.0,
                slat_width=20,
                ap_width=5,
                slant=75,
                n_inside=1.7)
'''


model = RayTraceModel(sources=[source],
                      optics=[test])
model.configure_traits()


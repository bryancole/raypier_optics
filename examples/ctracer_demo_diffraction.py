

import sys
sys.path.append('..')
import pyximport
pyximport.install()

from raytrace.tracer import RayTraceModel
from raytrace.sources import BroadbandRaySource
from raytrace.diffraction_gratings import RectangularGrating


source = BroadbandRaySource(number=5,
                            wavelength_start=1.5,
                            wavelength_end=1.6,
                            uniform_deltaf=True,
                            origin=(-50,0,0),
                            direction=(1.,0,0))

grating = RectangularGrating(length=25.4,
                             width=25.4,
                             thickness=6.0,
                             centre=(0.0,0.0,0.0),
                             direction=(-1,0,0),
                             lines_per_mm=600,
                             efficiency=1.0,
                             order=1
                             )

rays = source.InputRays

model = RayTraceModel(optics=[grating],
                    sources=[source,])

model.configure_traits()
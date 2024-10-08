

import sys
sys.path.append('..')
import pyximport
pyximport.install()

from raypier.tracer import RayTraceModel
from raypier.sources import BroadbandRaySource
from raypier.diffraction_gratings import RectangularGrating


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

rays = source.input_rays

model = RayTraceModel(optics=[grating],
                    sources=[source,])

model.configure_traits()
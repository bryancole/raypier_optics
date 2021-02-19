
from raypier.api import RayTraceModel, GaussletCapturePlane, EFieldPlane
from raypier.gausslet_sources import GaussianPointSource
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface

src = GaussianPointSource(origin=(0,0,0),
                          direction=(0,0,1),
                          E_vector=(1,0,0),
                          E1_amp = 1.0,
                          E2_amp = 0.0,
                          working_dist = 25.0,
                          numerical_aperture=0.1,
                          wavelength = 1.0,
                          beam_waist = 5.0,
                          display="wires",
                          opacity=0.1
                          )


cap = GaussletCapturePlane(centre=(0,0,30),
                           direction=(0,0,1),
                           width=10.,
                           height=10.)

p = EFieldPlane(detector=cap,
                align_detector=True,
                centre = (0,0,30),
                direction= (0,0,1),
                size=100,
                width=5.0,
                height=5.0)

img = IntensityImageView(field_probe=p)
surf = IntensitySurface(field_probe=p)

model = RayTraceModel(sources=[src], probes=[cap,p], results=[img,surf])
model.configure_traits()


 
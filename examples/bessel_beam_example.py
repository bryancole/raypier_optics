

from raypier.lenses import GeneralLens
from raypier.faces import AxiconFace, PlanarFace
from raypier.materials import OpticalMaterial
from raypier.shapes import CircleShape
from raypier.tracer import RayTraceModel
from raypier.gausslet_sources import CollimatedGaussletSource
from raypier.fields import EFieldPlane
from raypier.probes import GaussletCapturePlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface

shape = CircleShape(radius=2.0)

face1 = PlanarFace(z_height=0.0)
face2 = AxiconFace(z_height=1.0, gradient=0.1)

mat = OpticalMaterial(glass_name="N-BK7")

axicon = GeneralLens(name = "My Axicon",
                     centre = (0,0,0),
                     direction=(0,0,1),
                     shape=shape, 
                     surfaces=[face2, 
                               face1], 
                     materials=[mat])

src = CollimatedGaussletSource(origin=(0.001,0,-5.0),
                               direction=(0,0,1),
                               wavelength=0.5,
                               radius=1.0,
                               beam_waist=10.0,
                               resolution=10,
                               max_ray_len=50.0,
                               display='wires',
                               opacity=0.2
                               )

###Add some sensors
capture = GaussletCapturePlane(centre=(0,0,13), 
                               direction=(0,0,1),
                               width=5.0,
                               height=5.0)

field = EFieldPlane(centre=(0,0,13),
                    direction=(0,0,1),
                    detector=capture,
                    align_detector=True,
                    size=100,
                    width=2,
                    height=2)

image = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


model = RayTraceModel(optics=[axicon], 
                      sources=[src], 
                      probes=[capture,field],
                      results=[image,surf])
model.configure_traits()

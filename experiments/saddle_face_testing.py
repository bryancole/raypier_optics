

from raypier.tracer import RayTraceModel
from raypier.shapes import CircleShape
from raypier.faces import DistortionFace, PlanarFace, SphericalFace, SaddleFace
from raypier.general_optic import GeneralLens
from raypier.materials import OpticalMaterial
from raypier.distortions import SimpleTestZernikeJ7, NullDistortion
from raypier.gausslet_sources import CollimatedGaussletSource
from raypier.fields import EFieldPlane
from raypier.probes import GaussletCapturePlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface

shape = CircleShape(radius=10.0)
#f1 = SphericalFace(z_height=0.0, curvature=-25.0)
f1 = SaddleFace(z_height=0.0, curvature=0.0001)
f2 = PlanarFace(z_height=5.0)

#dist = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=0.01)
#dist = NullDistortion()
#f1 = DistortionFace(base_face=f1, distortion=dist)

mat = OpticalMaterial(glass_name="N-BK7")
lens = GeneralLens(direction=(0,0.0,1.0),
                   shape=shape, surfaces=[f1,f2], materials=[mat])

src = CollimatedGaussletSource(radius=8.0, resolution=1.4, 
                               origin=(0,0,-15), direction=(0,0,1),
                               wavelength=1.0,
                               beam_waist=15.0,
                               display="wires", opacity=0.2, show_normals=True)
src.max_ray_len=50.0


cap = GaussletCapturePlane(centre = (0,0,10),
                           direction= (0,0,1),
                           width=22,
                           height=22)

field = EFieldPlane(detector=cap,
                    align_detector=True,
                    centre=(0,0,10),
                    size=100,
                    width=20.0,
                    height=20.0)

img = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


model = RayTraceModel(optics=[lens], sources=[src], probes=[field,cap],
                      results=[img,surf])
model.configure_traits()


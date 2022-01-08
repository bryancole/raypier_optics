
from raypier.tracer import RayTraceModel
from raypier.faces import PlanarFace, ExtendedPolynomialFace, ConicFace
from raypier.shapes import CircleShape
from raypier.general_optic import GeneralLens
from raypier.materials import OpticalMaterial
from raypier.gausslet_sources import CollimatedGaussletSource
from raypier.fields import EFieldPlane
from raypier.probes import GaussletCapturePlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface

import numpy

coefs = numpy.array([[0.1, 0.1, 1],
                     [0.1, 0.2, 0.0],
                     [0.1, 0.0, 0.0]])

#coefs = numpy.zeros((3,3), 'd')


shape = CircleShape(radius=10.0)
f1 = ExtendedPolynomialFace(z_height=-2.0, curvature=-25.0, conic_const=0.0, norm_radius=5.0, coefs=coefs)
#f1 = ConicFace(z_height=0.0, curvature=-25.0, conic_const=0.0)
f2 = PlanarFace(z_height=5.0)

mat = OpticalMaterial(glass_name="N-BK7")
lens = GeneralLens(shape=shape, surfaces=[f1,f2], materials=[mat])

src = CollimatedGaussletSource(radius=8.0, resolution=6, 
                               origin=(0,0,-15), direction=(0,0,1),
                               display="wires", opacity=0.2, show_normals=True)
src.max_ray_len=50.0


cap = GaussletCapturePlane(centre = (0,0,50),
                           direction= (0,0,1),
                           width=20,
                           height=20)

field = EFieldPlane(detector=cap,
                    align_detector=True,
                    size=100,
                    width=1,
                    height=1)

img = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


model = RayTraceModel(optics=[lens], sources=[src], probes=[field,cap],
                      results=[img,surf])
model.configure_traits()










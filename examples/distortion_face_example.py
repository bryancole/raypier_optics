

from raypier.tracer import RayTraceModel

from raypier.shapes import CircleShape

from raypier.faces import DistortionFace, PlanarFace, SphericalFace

from raypier.general_optic import GeneralLens
from raypier.materials import OpticalMaterial
from raypier.distortions import SimpleTestZernikeJ7
from raypier.gausslet_sources import CollimatedGaussletSource


shape = CircleShape(radius=10.0)

#f1 = PlanarFace(z_height=0.0)
f1 = SphericalFace(z_height=0.0, curvature=-25.0)
f2 = PlanarFace(z_height=5.0)

dist = SimpleTestZernikeJ7(unit_radius=10.0, amplitude=1.0)

f1 = DistortionFace(base_face=f1, distortion=dist)

mat = OpticalMaterial(glass_name="N-BK7")

lens = GeneralLens(shape=shape, surfaces=[f1,f2], materials=[mat])

src = CollimatedGaussletSource(radius=3.0, resolution=1.5, 
                               origin=(0,0,-15), direction=(0,0,1),
                               display="wires", opacity=0.2, show_normals=True)
src.max_ray_len=50.0

model = RayTraceModel(optics=[lens], sources=[src])
model.configure_traits()


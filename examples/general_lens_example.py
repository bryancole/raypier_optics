



from raytrace.shapes import CircleShape, RectangleShape
from raytrace.faces import SphericalFace

from raytrace.materials import OpticalMaterial
from raytrace.lenses import GeneralLens
from raytrace.tracer import RayTraceModel


s1 = CircleShape(radius=10.0) ^ RectangleShape(width=5.,height=3.)

f1 = SphericalFace(z_height=4.0, curvature=50.0)
f2 = SphericalFace(z_height=-4.0, curvature=-50.0, invert=True)

lens = GeneralLens(shape=s1,
                   surfaces=[f1,f2])


model = RayTraceModel(optics=[lens])

model.configure_traits()




from raytrace.shapes import CircleShape, RectangleShape
from raytrace.faces import SphericalFace, PlanarFace

from raytrace.materials import OpticalMaterial
from raytrace.lenses import GeneralLens
from raytrace.tracer import RayTraceModel


#s1 = CircleShape(radius=10.0) ^ RectangleShape(width=5.,height=3.)
#s1 = RectangleShape(width=30,height=25) ^ CircleShape(radius=5.0)

s1 = RectangleShape(width=30,height=20) | RectangleShape(width=20,height=30) | CircleShape(radius=5.0, centre=(-10,-10)) |\
        CircleShape(radius=5.0, centre=(10,-10)) | CircleShape(radius=5.0, centre=(10,10)) | CircleShape(radius=5.0, centre=(-10,10))

f1 = SphericalFace(z_height=8.0, curvature=50.0)
f2 = PlanarFace(z_height=1.0)
f3 = SphericalFace(z_height=-8.0, curvature=-50.0, invert=True)
f4 = PlanarFace(z_height=-9.0, invert=True)

faces = [f1,f2,f3,f4]

lens = GeneralLens(shape=s1,
                   surfaces=faces)


model = RayTraceModel(optics=[lens])

model.configure_traits()
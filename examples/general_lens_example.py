



from raytrace.shapes import CircleShape, RectangleShape
from raytrace.faces import SphericalFace, PlanarFace

from raytrace.materials import OpticalMaterial
from raytrace.lenses import GeneralLens
from raytrace.tracer import RayTraceModel
from raytrace.sources import HexagonalRayFieldSource


#s1 = CircleShape(radius=20.0) ^ RectangleShape(width=5.,height=3.)
s1 = RectangleShape(width=30,height=25) ^ CircleShape(radius=5.0)

# s1 = RectangleShape(width=30,height=20) | RectangleShape(width=20,height=30) | CircleShape(radius=5.0, centre=(-10,-10)) |\
#         CircleShape(radius=5.0, centre=(10,-10)) | CircleShape(radius=5.0, centre=(10,10)) | CircleShape(radius=5.0, centre=(-10,10))

f1 = SphericalFace(z_height=8.0, curvature=50.0)
m1 = OpticalMaterial(from_database=False, refractive_index=1.5)
f2 = PlanarFace(z_height=1.0)
m2 = OpticalMaterial(from_database=False, refractive_index=1.1)
f3 = SphericalFace(z_height=-8.0, curvature=-50.0, invert=True)
m3 = OpticalMaterial(from_database=False, refractive_index=1.6)
f4 = PlanarFace(z_height=-9.0, invert=False)

faces = [f1,f2,f3,f4]
mats = [m1,m2,m3]

lens = GeneralLens(centre=(0,0,50),
                    shape=s1,
                   surfaces=faces,
                   materials=mats)

src = HexagonalRayFieldSource(gauss_width=5.0,
                              display="wires",
                              opacity=0.1,
                              show_normals=True)


model = RayTraceModel(optics=[lens], sources=[src])

model.configure_traits()




from raypier.api import CircleShape, RectangleShape, PolygonShape, HexagonShape
from raypier.api import SphericalFace, PlanarFace

from raypier.api import OpticalMaterial, OpticalMaterialFromFormula
from raypier.api import GeneralLens
from raypier.api import RayTraceModel
from raypier.api import HexagonalRayFieldSource


#s1 = RectangleShape(width=30,height=25) ^ CircleShape(radius=3.0)
# s1 = PolygonShape(coordinates=list(reversed([(-30.0,0.),
#                                (-15.,25.98),
#                                (15.0,25.98),
#                                (30.,0.),
#                                (15.,-25.98),
#                                (-15.,-25.98)])))
s1 = HexagonShape(radius=15)
s1 = s1 ^ CircleShape(radius=3.0)

f1 = SphericalFace(z_height=8.0, curvature=50.0)
m1 = OpticalMaterial(from_database=False, refractive_index=1.5, name="Mat1")
f2 = PlanarFace(z_height=1.0)
m2 = OpticalMaterial(from_database=False, refractive_index=1.1, name="Mat2")
f3 = SphericalFace(z_height=-8.0, curvature=-50.0)
m3 = OpticalMaterial(from_database=False, refractive_index=1.6, name="Mat3")
f4 = PlanarFace(z_height=-9.0, invert=False)

###If you want to specify a custom material you can do this:
# m3 = OpticalMaterialFromFormula(formula_id=1,
#                                 coefs=[0.0,
#                                 0.6961663,
#                                 0.0684043,
#                                 0.4079426,
#                                 0.1162414,
#                                 0.8974794,
#                                 9.896161])

faces = [f1,f2,f3,f4]
mats = [m1, m2,m3]

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
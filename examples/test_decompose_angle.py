

from raytrace.gausslet_sources import SingleGaussletSource, CollimatedGaussletSource
from raytrace.gausslets import AngleDecompositionPlane
from raytrace.lenses import PlanoConvexLens, GeneralLens
from raytrace.achromats import EdmundOptic45806
from raytrace.shapes import CircleShape
from raytrace.faces import SphericalFace
from raytrace.materials import OpticalMaterial
from raytrace.fields import EFieldPlane
from raytrace.intensity_image import IntensityImageView
from raytrace.intensity_surface import IntensitySurface
from raytrace.probes import GaussletCapturePlane

from raytrace.tracer import RayTraceModel

lens3 = EdmundOptic45806(centre=(0,0,50),
                         direction=(0,0,1))

lens2 = PlanoConvexLens(centre=(0,0,50),
                       n_inside=1.5,
                       curvature=50.0,
                       diameter=30.0)

shape = CircleShape(radius=12.3)
f1 = SphericalFace(curvature=50.0, z_height=2.5)
f2 = SphericalFace(curvature=-50.0, z_height=-2.5)
m = OpticalMaterial(glass_name="N-BK7")
lens = GeneralLens(centre=(0,0,50),
                    shape=shape, 
                    surfaces=[f1,f2], 
                    materials=[m],
                    coating_material=OpticalMaterial(from_database=False,refractive_index=1.25),
                    coating_thickness=0.2)



#src = SingleGaussletSource(beam_waist=1.0, max_ray_len=55.0)
src = CollimatedGaussletSource(radius=1.0, resolution=10,
                               beam_waist=10.0,
                               wavelength=1.0,
                               display="wires",
                               opacity=0.05)

cap = GaussletCapturePlane(centre=(0,0,120))

probe = EFieldPlane(detector=cap,
                    align_detector=True,
                    source=src,
                    centre=(0,0,130),
                    direction=(0,0,1),
                    exit_pupil_offset=0.,
                    width=15,
                    height=15,
                    size=100)

img = IntensityImageView(field_probe=probe)
surf = IntensitySurface(field_probe=probe)

decomp = AngleDecompositionPlane(sample_spacing=5.0,
                                 width=128, height=128,
                                 centre = (0,0,100.52),
                                 direction = (0,0,1),
                                 diameter=25.0,
                                 max_angle=2.0)


model = RayTraceModel(sources=[src], optics=[lens2, decomp], probes=[probe, cap],
                      results=[img, surf])
model.configure_traits()


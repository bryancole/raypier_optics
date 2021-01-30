
import faulthandler
faulthandler.enable()


from raypier.gausslet_sources import SingleGaussletSource, CollimatedGaussletSource
from raypier.lenses import PlanoConvexLens, GeneralLens
from raypier.shapes import CircleShape
from raypier.faces import SphericalFace
from raypier.materials import OpticalMaterial
from raypier.fields import EFieldPlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface
from raypier.probes import GaussletCapturePlane

from raypier.tracer import RayTraceModel

lens = PlanoConvexLens(centre=(0,0,50),
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
                    materials=[m])



#src = SingleGaussletSource(beam_waist=1.0, max_ray_len=55.0)
src = CollimatedGaussletSource(radius=10.0, resolution=10,
                               beam_waist=10.0,
                               wavelength=1.0)

cap = GaussletCapturePlane(centre=(0,0,120))

probe = EFieldPlane(detector=cap,
                    align_detector=True,
                    source=src,
                    centre=(0,0,70),
                    direction=(0,0,1),
                    exit_pupil_offset=0.,
                    width=15,
                    height=15,
                    size=100)

img = IntensityImageView(field_probe=probe)
surf = IntensitySurface(field_probe=probe)



model = RayTraceModel(sources=[src], optics=[lens], probes=[probe, cap],
                      results=[img, surf])
model.configure_traits()
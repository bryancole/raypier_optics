

from raypier.tracer import RayTraceModel
from raypier.shapes import CircleShape
from raypier.faces import DistortionFace, PlanarFace, SphericalFace
from raypier.general_optic import GeneralLens
from raypier.materials import OpticalMaterial
from raypier.distortions import SimpleTestZernikeJ7, NullDistortion, ZernikeSeries
from raypier.gausslet_sources import CollimatedGaussletSource
from raypier.fields import EFieldPlane
from raypier.probes import GaussletCapturePlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface
from raypier.api import Constraint

from traits.api import Range, observe
from traitsui.api import View, Item, VGroup


shape = CircleShape(radius=10.0)

#f1 = SphericalFace(z_height=0.0, curvature=-25.0)

f1 = PlanarFace(z_height=0.0)
f2 = PlanarFace(z_height=5.0)

dist = ZernikeSeries(unit_radius=10.0, coefficients=[(i,0.0) for i in range(12)])
f1 = DistortionFace(base_face=f1, distortion=dist)


class Sliders(Constraint):
    """Make a Constrain object just to give us a more convenient UI for adjusting Zernike coefficients.
    """
    J0 = Range(-1.0,1.0,0.0)
    J1 = Range(-1.0,1.0,0.0)
    J2 = Range(-1.0,1.0,0.0)
    J3 = Range(-1.0,1.0,0.0)
    J4 = Range(-1.0,1.0,0.0)
    J5 = Range(-1.0,1.0,0.0)
    J6 = Range(-1.0,1.0,0.0)
    J7 = Range(-1.0,1.0,0.0)
    J8 = Range(-1.0,1.0,0.0)
    
    traits_view = View(VGroup(
            Item("J0", style="custom"),
            Item("J1", style="custom"),
            Item("J2", style="custom"),
            Item("J3", style="custom"),
            Item("J4", style="custom"),
            Item("J5", style="custom"),
            Item("J6", style="custom"),
            Item("J7", style="custom"),
            Item("J8", style="custom"),
        ),
        resizable=True)
    
    def _anytrait_changed(self):
        dist.coefficients = list(enumerate([self.J0, self.J1, self.J2, self.J3, self.J4, 
                                            self.J5, self.J6, self.J7, self.J8]))


mat = OpticalMaterial(glass_name="N-BK7")
lens = GeneralLens(shape=shape, surfaces=[f1,f2], materials=[mat])

src = CollimatedGaussletSource(radius=9.0, resolution=20, 
                               origin=(0,0,-15), direction=(0,0,1),
                               display="wires", opacity=0.02,
                               wavelength=1.0, 
                               beam_waist=10.0,
                               show_normals=True)
src.max_ray_len=50.0


cap = GaussletCapturePlane(centre = (0,0,50),
                           direction= (0,0,1),
                           width=20,
                           height=20)

field = EFieldPlane(centre=(0,0,30),
                    detector=cap,
                    align_detector=True,
                    size=100,
                    width=20,
                    height=20)

img = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


model = RayTraceModel(optics=[lens], sources=[src], probes=[field,cap],
                      results=[img,surf], constraints=[Sliders()])
model.configure_traits()


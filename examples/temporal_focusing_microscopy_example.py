
from pathlib import Path

from traits.api import Range, Float, observe
from traitsui.api import View, Item, VGroup

from raypier.beamsplitters import UnpolarisingBeamsplitterCube
from raypier.lenses import PlanoConvexLens
from raypier.diffraction_gratings import RectangularGrating
from raypier.gausslet_sources import BroadbandGaussletSource, SingleGaussletSource
from raypier.tracer import RayTraceModel
from raypier.fields import EFieldPlane
from raypier.intensity_image import IntensityImageView
from raypier.intensity_surface import IntensitySurface
from raypier.probes import GaussletCapturePlane
from raypier.constraints import Constraint
from raypier.editors import NumEditor
from raypier.shapes import CircleShape
from raypier.general_optic import GeneralLens
from raypier.faces import AsphericFace, SphericalFace
from raypier.materials import OpticalMaterial



src = BroadbandGaussletSource(
    origin = (0,0,0),
    direction=(1.0,0.0,0.0),
    E_vector=(0,0,1),
    working_dist=0.0,
    number=200,
    wavelength = 1.0,
    wavelength_extent = 0.03,
    bandwidth_nm = 13.0,
    beam_waist = 1000.0,
    display='wires',
    show_paras=False
    )

# src = SingleGaussletSource(
#     origin = (0,0,0),
#     direction=(1.,0.,0.0),
#     working_dist=0.0,
#     wavelength=0.8,
#     beam_waist=1.0, #in microns
#     E_vector=(0,0,1)
#     )

grating = RectangularGrating(centre=(220.0,0.,0.),
                             direction=(-1,1,0.),
                             length=15.,
                             width=20.0,
                             thickness=3.0,
                             lines_per_mm=1400.0)
grating.orientation = 45.5

lens1 = PlanoConvexLens(centre=(40.0,0.0,0.0),
                        direction=(-1,0,0),
                        diameter=25.0,
                        CT=6.0,
                        n_inside=1.6,
                        curvature=100.0)

### Not used now. Use Aspheric instead.
lens2 = PlanoConvexLens(centre=(10.0,-30.0,0.0),
                        direction=(0,1,0),
                        diameter=25.0,
                        CT=6.0,
                        n_inside=1.6,
                        curvature=25.0)

### Construct an aspheric objective using General Lens framework.
### This is intended to represent Edmund Optics #66-330
circle = CircleShape(radius=10.0)
s1 = AsphericFace(z_height=0.0,
                  curvature=1./8.107287E-02,
                  conic_const=-6.196140E-01,
                  A6=-1.292772E-08, #I'm not 100% sure if the signs of the polynomial terms are right.
                  A8=-1.932447E-10
                  )
s2 = SphericalFace(curvature=-200.0,
                   z_height=-8.0)
mat = OpticalMaterial(glass_name="L-BAL35")
asphere = GeneralLens(name="Sample Objective",
                      centre=(10.0,-30.0,0.0),
                      direction=(0,1,0),
                      shape=circle,
                      surfaces=[s2,s1],
                      materials=[mat])

bs = UnpolarisingBeamsplitterCube(centre = (10.0, 0., 0.),
                                  size=15.0,
                                  )

capture = GaussletCapturePlane(centre=(10,-53.3,0), 
                               direction=(0,1,0),
                               width=15.0,
                               height=15.0)

field = EFieldPlane(centre=(10,-53.3,0),
                    direction=(0,0,1),
                    detector=capture,
                    align_detector=True,
                    size=100,
                    width=0.1,
                    height=0.5,
                    time_ps=-7.0)

image = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


class MyConstraint(Constraint):
    time = Range(-15.0,15.0,-7.0)
    time_offset = Float(0.0)
    
    traits_view = View(VGroup(
        Item("time", style='custom'),
        Item("time_offset", editor=NumEditor)
        ))
    
    @observe("time, time_offset")
    def on_time_change(self, evt):
        field.time_ps = self.time + self.time_offset
        
        
    def animate(self, dt, count):
        for i in range(count):
            field.time_ps = self.time + self.time_offset + i*dt
            U = field.intensity
            image.save_plot(f"{Path.home()}/range_{i:03d}.png")

cst = MyConstraint()

model = RayTraceModel(optics = [bs, grating, lens1, asphere],
                      sources = [src], probes=[field, capture],
                      results=[image,surf],
                      constraints = [cst])

model.configure_traits()


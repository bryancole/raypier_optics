
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


src = BroadbandGaussletSource(
    origin = (0,0,0),
    direction=(1,0,0.001),
    E_vector=(0,0,1),
    working_dist=0.0,
    wavelength = 1.0,
    wavelength_extent = 0.01,
    bandwidth_nm = 1.0,
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

lens2 = PlanoConvexLens(centre=(10.0,-30.0,0.0),
                        direction=(0,1,0),
                        diameter=25.0,
                        CT=6.0,
                        n_inside=1.6,
                        curvature=25.0)

bs = UnpolarisingBeamsplitterCube(centre = (10.0, 0., 0.),
                                  size=15.0,
                                  )

capture = GaussletCapturePlane(centre=(10,-66,0), 
                               direction=(0,1,0),
                               width=15.0,
                               height=15.0)

field = EFieldPlane(centre=(10,-66,0),
                    direction=(0,0,1),
                    detector=capture,
                    align_detector=True,
                    size=100,
                    width=5,
                    height=5)

image = IntensityImageView(field_probe=field)
surf = IntensitySurface(field_probe=field)


class MyConstraint(Constraint):
    time = Range(-15.0,15.0,0.0)
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
            print("Calc:", i)
            image.save_plot(f"/home/bryan/tmp/range_{i:03d}.png")

cst = MyConstraint()

model = RayTraceModel(optics = [bs, grating, lens1, lens2],
                      sources = [src], probes=[field, capture],
                      results=[image,surf],
                      constraints = [cst])

model.configure_traits()


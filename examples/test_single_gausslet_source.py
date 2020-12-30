
from raytrace.tracer import RayTraceModel
from raytrace.sources import SingleGaussletSource
from raytrace.lenses import PlanoConvexLens
from raytrace.apertures import CircularAperture
from raytrace.fields import EFieldPlane
from raytrace.constraints import BaseConstraint
from raytrace.intensity_image import IntensityImageView
from raytrace.intensity_surface import IntensitySurface

from traits.api import Range, on_trait_change
from traitsui.api import View, Item


src = SingleGaussletSource(direction=(0,0,1),
                              radius=1000.0,
                              wavelength=100.0)

src.InputRays

probe = EFieldPlane(source=src,
                    centre=(0,0,50),
                    direction=(0,0,1),
                    exit_pupil_offset=0.,
                    width=5.0,
                    height=5.0,
                    size=100)

img = IntensityImageView(field_probe=probe)
surf = IntensitySurface(field_probe=probe)


class FocalPlane(BaseConstraint):
    z_pos = Range(50.0,130.0, 57.73)
    
    traits_view = View(Item("z_pos"))
    
    def __init__(self, *args, **kwds):
        super().__init__(*args, **kwds)
        self.on_change_z_pos()
    
    @on_trait_change("z_pos")
    def on_change_z_pos(self):
        probe.centre = (0,0,self.z_pos)
        
        

model = RayTraceModel(sources=[src], optics=[],
                      probes=[probe], constraints=[FocalPlane()],
                      results=[img, surf])



model.configure_traits()


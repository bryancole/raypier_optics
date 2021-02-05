"""
Traits wrappers for Distortion objects
"""

from traits.api import Instance, Float, HasTraits, observe, Event
from traitsui.api import View, VGroup, Item

from .editors import NumEditor
from .core import cdistortions, ctracer


class BaseDistortion(HasTraits):
    updated = Event()
    c_distortion = Instance(ctracer.Distortion, ())


class SimpleTestZernikeJ7(BaseDistortion):
    unit_radius = Float(10.0)
    amplitude = Float(0.001)
    
    c_distortion = Instance(cdistortions.SimpleTestZernikeJ7, ())
    
    traits_view = View(VGroup(
                    Item("unit_radius", editor=NumEditor),
                    Item("amplitude", editor=NumEditor)
                        )
                    )
    
    def _c_distortion_default(self):
        return cdistortions.SimpleTestZernikeJ7(unit_radius=self.unit_radius,
                                                amplitude=self.amplitude)
        
    @observe("unit_radius, amplitude")
    def on_params_change(self, evt):
        d = self.c_distortion
        setattr(d, evt.name, evt.new)
        self.updated = True
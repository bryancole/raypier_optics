"""
Traits wrappers for Distortion objects
"""

from traits.api import Instance, Float, HasTraits, observe, Event, List, Tuple, Int
from traitsui.api import View, VGroup, Item
from traitsui.editors import ListEditor, TupleEditor

from .editors import NumEditor
from .core import cdistortions, ctracer
from traits.observation.events import ListChangeEvent
from traits.observation._trait_change_event import TraitChangeEvent


class BaseDistortion(HasTraits):
    updated = Event()
    c_distortion = Instance(ctracer.Distortion, ())
    
    
class NullDistortion(BaseDistortion):
    pass


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
        
        
coef_editor = TupleEditor(cols=2, labels=["j", "value"])
        
        
class ZernikeSeries(BaseDistortion):
    unit_radius = Float(10.0)
    
    coefficients = List(Tuple(Int,Float))
    
    c_distortion = Instance(cdistortions.ZernikeDistortion)
    
    traits_view = View(VGroup(
                    Item("unit_radius", editor=NumEditor),
                    Item("coefficients", editor=ListEditor(editor=coef_editor))
                        )
                    )
    
    def __init__(self, *args, **kwds):
        coefs = dict(kwds.get("coefficients", []))
        if args:
            coefs.update(args[0])
        jks = {k for k in kwds if k.startswith("j")}
        coefs.update({int(k[1:]):float(kwds[k]) for k in jks})
        for k in jks:
            del kwds[k]
        kwds['coefficients'] = sorted(coefs.items())
        super().__init__(**kwds)
    
    def _c_distortion_default(self):
        return cdistortions.ZernikeDistortion(self.coefficients,
                                              unit_radius=self.unit_radius)
        
    @observe("coefficients.items, unit_radius")
    def on_params_change(self, evt):
        if isinstance(evt, TraitChangeEvent) and evt.name == "unit_radius":
            self.c_distortion.unit_radius = evt.new
        else:
            self.c_distortion.set_coefs(self.coefficients)
        self.updated = True
        
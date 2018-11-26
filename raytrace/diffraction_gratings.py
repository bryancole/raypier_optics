'''
Created on 26 Nov 2018

@author: bryan
'''
from traits.api import Float, Int, on_trait_change

from traitsui.api import View, Item, VGroup

from raytrace.cmaterials import DiffractionGratingMaterial

from raytrace.bases import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor, Optic
from raytrace.mirrors import RectMirror


class RectangularGrating(RectMirror):
    name = "Rectangular Diffraction Grating"
    order = Int(1)
    lines_per_mm = Float(600.0)
    efficiency = Float(1.0)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('length', editor=NumEditor),
                       Item('width', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('offset', editor=NumEditor),
                       Item('order'),
                       Item('lines_per_mm', editor=NumEditor),
                       Item('efficiency', editor=NumEditor)
                        ),
                   )
    
    def _material_default(self):
        return DiffractionGratingMaterial(order=self.order,
                                          lines_per_mm=self.lines_per_mm,
                                          efficiency=self.efficiency)
        
    @on_trait_change("order, lines_per_mm, efficiency")
    def on_parameter_change(self, name, new):
        setattr(self.material, name, new)
        self.update = True
        
        
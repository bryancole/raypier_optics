'''
Created on 26 Nov 2018

@author: bryan
'''
from traits.api import Float, Int, on_trait_change

from traitsui.api import View, Item, VGroup

from raypier.core.cmaterials import DiffractionGratingMaterial

from raypier.bases import Traceable, normaliseVector, NumEditor,\
     ComplexEditor, VectorEditor, Optic
from raypier.mirrors import RectMirror


class RectangularGrating(RectMirror):
    """
    A reflective diffraction-grating.
    """
    
    #: Object name is display in model-tree
    name = "Rectangular Diffraction Grating"
    
    #: Order of diffraction (int).
    order = Int(1)
    
    #: Lines per mm for the grating
    lines_per_mm = Float(600.0)
    
    #: The diffraction efficiency. A float in range 0.0 -> 1.0 (no loss).
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
    
    @on_trait_change("centre")
    def on_centre_changed(self, vnew):
        self.material.origin = vnew
    
    def _material_default(self):
        return DiffractionGratingMaterial(order=self.order,
                                          origin=self.centre,
                                          lines_per_mm=self.lines_per_mm,
                                          efficiency=self.efficiency)
        
    @on_trait_change("order, lines_per_mm, efficiency")
    def on_parameter_change(self, name, new):
        setattr(self.material, name, new)
        self.update = True
        
    def make_step_shape(self):
        shape, colour = super().make_step_shape()
        return shape, "pink"
        
        
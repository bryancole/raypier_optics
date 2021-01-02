
"""
Traited wrappers for cfaces objects
"""


from .core.cfaces import ShapedPlanarFace, ShapedSphericalFace, ConicRevolutionFace, AsphericFace as AsphericFace_
from .shapes import BaseShape
from .editors import NumEditor

from traits.api import HasStrictTraits, Instance, Float, Bool, observe
from traitsui.api import View, VGroup, Item, VGrid


class PlanarFace(HasStrictTraits):
    z_height = Float(0.0)
     
    cface = Instance(ShapedPlanarFace, (), )
    
    traits_view = View(VGroup(Item("z_height", editor=NumEditor, tooltip="surface height in mm")))
    
    def _z_height_changed(self, znew):
        self.cface.z_plane = znew
        
        
class SphericalFace(PlanarFace):
    curvature = Float(100.0)
    
    cface = Instance(ShapedSphericalFace, ())
    
    traits_view = View(VGroup(
            Item("z_height", editor=NumEditor, tooltip="surface height at x=0,y=0 in mm"),
            Item("curvature", editor=NumEditor, tooltip="radius of curvature, in mm")
        ))
    
    def _curvature_changed(self, cnew):
        self.cface.curvature = cnew
        
        
class ConicFace(SphericalFace):
    conic_const = Float(0.0)
    invert = Bool(False)
    
    cface = Instance(ConicRevolutionFace, ())
    
    traits_view = View(VGroup(
            Item("z_height", editor=NumEditor, tooltip="surface height at x=0,y=0 in mm"),
            Item("curvature", editor=NumEditor, tooltip="radius of curvature, in mm"),
            Item("conic_const", editor=NumEditor, tooltip="conic constant for the surface"),
            Item("invert", editor=NumEditor, tooltip="Invert the outward sense for this surface")
        ))
    
    def _conic_const_changed(self, knew):
        self.cface.conic_const = knew
    
    def _invert_changed(self, vnew):
        self.cface.invert_normals = int(vnew)
    
    
class AsphericFace(ConicFace):
    A4 = Float(0.0)
    A6 = Float(0.0)
    A8 = Float(0.0)
    A10 = Float(0.0)
    A12 = Float(0.0)
    A14 = Float(0.0)
    A16 = Float(0.0)
                
    cface = Instance(AsphericFace_, ())
    
    traits_view = View(VGroup(
            Item("z_height", editor=NumEditor, tooltip="surface height at x=0,y=0 in mm"),
            Item("curvature", editor=NumEditor, tooltip="radius of curvature, in mm"),
            Item("conic_const", editor=NumEditor, tooltip="conic constant for the surface"),
            VGrid(
                Item("A4", editor=NumEditor), Item("A6", editor=NumEditor), Item("A8", editor=NumEditor), 
                Item("A10", editor=NumEditor), Item("A12", editor=NumEditor),
                Item("A14", editor=NumEditor), Item("A16", editor=NumEditor),
                show_border=True,
                label="Polynomial coefs"
                ),
            Item("invert", tooltip="Invert the outward sense for this surface")
        ))
    
    @observe("A4, A6, A8, A10, A12, A14, A16")
    def on_coef_change(self, evt):
        setattr(self.cface, evt.name, evt.new)
        
    
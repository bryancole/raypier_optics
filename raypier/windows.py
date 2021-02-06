

from traits.api import Float, Str, Instance
from traitsui.api import Item, View, VGroup

from tvtk.api import tvtk

from .core.cfaces import CircularFace
from .core.cmaterials import CoatedDispersiveMaterial
from .dispersion import NamedDispersionCurve, NondispersiveCurve
from .bases import Optic, Traceable, NumEditor
from .core.ctracer import FaceList


class CircularWindow(Optic):
    abstract = False
    
    thickness = Float(5.0, desc="centre thickness")
    diameter = Float(15.0)
    offset = Float(0.0)
    
    glass_type = Str("N-BK7")
    coating_thickness = Float(0.283) #microns
    coating_refractive_index = Float(1.37)
    
    vtk_cylinder = Instance(tvtk.Cylinder, ())
    
    traits_view = View(VGroup(
                   Traceable.uigroup,  
                   Item('glass_type'),
                   Item('thickness', editor=NumEditor),
                   Item('diameter', editor=NumEditor),
                   Item('offset', editor=NumEditor)
                   )
                )
    

    def _material_default(self):
        glass = NamedDispersionCurve(self.glass_type)
        cri = self.coating_refractive_index
        coating = NondispersiveCurve(refractive_index=cri)
        air = NondispersiveCurve(1.0)
        
        m = CoatedDispersiveMaterial(dispersion_inside=glass,
                                     dispersion_outside=air,
                                     dispersion_coating=coating)
        return m
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        fl.faces = [CircularFace(owner=self, diameter=self.diameter,
                                 offset = self.offset, z_plane=0.0,
                                material = self.material), 
                    CircularFace(owner=self, diameter=self.diameter,
                                material=self.material, offset=self.offset,
                                z_plane=self.thickness)]
        return fl
    
    

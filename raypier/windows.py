

from traits.api import Float, Str, Instance, on_trait_change
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
    
    vtk_cylinder = Instance(tvtk.CylinderSource, (), {'resolution':32})
    cyl_trans = Instance(tvtk.Transform, (), transient=True)
    
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
                                z_plane=self.thickness, invert_normals=True)]
        return fl
    
    @on_trait_change("diameter, thickness, offset")
    def config_pipeline(self):
        cyl = self.vtk_cylinder
        cyl.radius = self.diameter/2.
        thick = self.thickness
        cyl.height = thick
        
        t = self.cyl_trans
        t.identity()
        t.translate(self.offset,0,thick/2.)
        t.rotate_x(90.)
        
        self.update = True
    
    def _pipeline_default(self):
        cyl = self.vtk_cylinder
        norms = tvtk.PolyDataNormals(input_connection=cyl.output_port)
        transF1 = tvtk.TransformFilter(input_connection=norms.output_port, 
                                       transform=self.cyl_trans)
        transF2 = tvtk.TransformFilter(input_connection=transF1.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return transF2
    

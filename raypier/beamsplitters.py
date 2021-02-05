
from traits.api import Float, on_trait_change, observe

from traitsui.api import View, VGroup, Item

from .bases import Traceable
from .prisms import Extrusion
from .editors import NumEditor, ComplexEditor
from .core.cfaces import ExtrudedPlanarFace
from .core.ctracer import FaceList
from .core.cmaterials import FullDielectricMaterial, \
            LinearPolarisingMaterial, PartiallyReflectiveMaterial, FullDielectricMaterial,\
            SingleLayerCoatedMaterial
            
import numpy



class BaseBeamsplitterCube(Extrusion):
    abstract = True
    size = Float(10.0)
    name = "BaseBeamsplitterCube"
    
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('size', editor=NumEditor),
                       Item('n_inside', editor=ComplexEditor),
                       )
                       )

    @on_trait_change("size")
    def config_profile(self):
        h = self.size/2
        self.z_height_1 = -h
        self.z_height_2 = h
        self.profile = [(-h,-h), (-h,h), (h,h), (h,-h)]


class PolarisingBeamsplitterCube(BaseBeamsplitterCube):
    abstract=False
    name = "Polarising Cube Beamsplitter"
    
    def make_faces(self):
        h = self.size/2
        faces = super(PolarisingBeamsplitterCube,self).make_faces()
        m = LinearPolarisingMaterial()
        p_face = ExtrudedPlanarFace(owner=self, 
                                    z1=self.z_height_1, 
                                    z2=self.z_height_2, 
                                    x1=-h, y1=-h, 
                                    x2=h, y2=h, 
                                    material=m)
        faces.append( p_face )
        return faces
    
    def _vtk_profile(self):
        h = self.size/2
        profile = numpy.concatenate([ self.profile, [(-h,-h), (h,h)] ])
        return profile
    
    
class UnpolarisingBeamsplitterCube(PolarisingBeamsplitterCube):
    abstract=False
    name = "Non_Polarising_Cube_Beamsplitter"
    
    reflectivity = Float(0.5)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('size', editor=NumEditor),
                       Item('reflectivity', editor=NumEditor),
                       Item('n_inside', editor=ComplexEditor),
                       )
                       )
    
    def _material_default(self):
        #print("n_inside", self.n_inside, "n_outside", self.n_outside)
        m = FullDielectricMaterial(n_inside = self.n_inside,
                                    n_outside = self.n_outside)
#         m = SingleLayerCoatedMaterial(n_inside = self.n_inside,
#                                     n_outside = self.n_outside)
        #print("thresholds:", m.reflection_threshold, m.transmission_threshold)
        return m
    
    def make_faces(self):
        h = self.size/2
        faces = super(PolarisingBeamsplitterCube,self).make_faces()
        m = PartiallyReflectiveMaterial(reflectivity=self.reflectivity)
        p_face = ExtrudedPlanarFace(owner=self, 
                                    z1=self.z_height_1, 
                                    z2=self.z_height_2, 
                                    x1=-h, y1=-h, 
                                    x2=h, y2=h, 
                                    material=m)
        faces.append( p_face )
        return faces
    
    @observe("reflectivity")
    def on_refl_changed(self, evt):
        refl = self.reflectivity
        self.faces.faces[-1].material.reflectivity = refl
        
    
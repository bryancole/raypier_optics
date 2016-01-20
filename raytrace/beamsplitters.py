
from traits.api import Float, on_trait_change

from traitsui.api import View, VGroup, Item

from raytrace.bases import Traceable
from raytrace.prisms import Extrusion
from raytrace.cfaces import ExtrudedPlanarFace
from raytrace.ctracer import FaceList
from raytrace.cmaterials import FullDielectricMaterial, \
            LinearPolarisingMaterial


class BaseBeamsplitterCube(Extrusion):
    abstract = True
    size = Float(10.0)
    name = "BaseBeamsplitterCube"
    
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('size'),
                       Item('n_inside'),
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
    
    
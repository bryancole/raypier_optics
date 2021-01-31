#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raypier.
#
#    Raypier is free software: you can redistribute it and/or modify
#    it under the terms of the GNU General Public License as published by
#    the Free Software Foundation, either version 3 of the License, or
#    (at your option) any later version.
#
#    This program is distributed in the hope that it will be useful,
#    but WITHOUT ANY WARRANTY; without even the implied warranty of
#    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
#    GNU General Public License for more details.
#
#    You should have received a copy of the GNU General Public License
#    along with this program.  If not, see <http://www.gnu.org/licenses/>.


from traits.api import HasTraits, Array, Float, Complex,\
            Property, List, Instance, on_trait_change, Range, Any,\
            Tuple, Event, cached_property, Set, Int, Trait, Bool
from traitsui.api import View, Item, ListEditor, VSplit,\
            RangeEditor, ScrubberEditor, HSplit, VGroup, Heading
from tvtk.api import tvtk
from tvtk.pyface.scene_model import SceneModel
from tvtk.pyface.scene_editor import SceneEditor
import numpy
from itertools import chain, islice, tee

#from raypier.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
#             Traceable, NumEditor, dotprod, transformPoints, transformNormals

from raypier.bases import Optic, Traceable, NumEditor, RaypierObject
from raypier.core.cfaces import ElipticalPlaneFace, CircularFace
from raypier.core.ctracer import FaceList
from raypier.mirrors import BaseMirror
from raypier.core.cmaterials import PECMaterial, FullDielectricMaterial


class HollowRetroreflector(BaseMirror):
    abstract = False
    name = "Hollow Retroreflector"
    diameter = Float(25.4)

    imp_cyl = Instance(tvtk.Cylinder, (), transient=True)    
    scale = Instance(tvtk.Transform, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                        ),
                   )
    
    def _faces_default(self):
        fl = FaceList(owner=self)
        sqrt2 = numpy.sqrt(2)
        x = sqrt2*numpy.cos(numpy.pi*2./3)
        y = sqrt2*numpy.sin(numpy.pi*2./3)
        m = self.material
        fl.faces = [ElipticalPlaneFace(owner=self, g_x=sqrt2, g_y=0, material=m),
                    ElipticalPlaneFace(owner=self, g_x=x, g_y=y, material=m),
                    ElipticalPlaneFace(owner=self, g_x=x, g_y=-y, material=m),
                    ]
        return fl
    
    @on_trait_change("diameter")
    def config_pipeline(self):
        r = self.diameter/2.
        cyl = self.imp_cyl
        cyl.radius = r
        s=self.scale
        s.identity()
        s.scale(r,r,r)
        
        self.update = True
    
    def _pipeline_default(self):
        sqrt2 = numpy.sqrt(2)
        x = -sqrt2*numpy.cos(numpy.pi*2./3)
        y = sqrt2*numpy.sin(numpy.pi*2./3)
        p1 = tvtk.PlaneSource(origin=(0,0,0), \
                              point1=(-sqrt2,0,1), point2=(x,y,1))
        p2 = tvtk.PlaneSource(origin=(0,0,0), \
                               point1=(x,y,1), point2=(x,-y,1))
        p3 = tvtk.PlaneSource(origin=(0,0,0), \
                               point1=(x,-y,1), point2=(-sqrt2,0,1))
        for p in (p1,p2,p3):
            p.set_resolution(32,32)
        append = tvtk.AppendPolyData()
        append.add_input_connection(p1.output_port)
        append.add_input_connection(p2.output_port)
        append.add_input_connection(p3.output_port)
        
        scale_f = tvtk.TransformFilter(input_connection=append.output_port, 
                                       transform=self.scale)
        
        trans = tvtk.Transform()
        trans.rotate_x(90.0)
        self.imp_cyl.transform = trans
        
        clip = tvtk.ClipPolyData(input_connection=scale_f.output_port, 
                                 clip_function = self.imp_cyl,
                                 inside_out=True)
                               
        transF2 = tvtk.TransformFilter(input_connection=clip.output_port, 
                                       transform=self.transform)
        self.config_pipeline()
        return transF2
    
    
class SolidRetroreflector(Optic, HollowRetroreflector):
    abstract = False
    yaml_tag = "!SolidRetroreflector"
    n_inside = 1.5
    name = "Solid Retroreflector"
    
    offset = 0.0
    
    thickness = Float(20.0, desc="distance from input face to apex")
    
    cyl = Instance(tvtk.CylinderSource, (), transient=True)
    cube = Instance(tvtk.CubeSource, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('n_inside'),
                        ),
                   )
    
    def _material_default(self):
        print("Making full dielectric material")
        m = FullDielectricMaterial(n_inside = self.n_inside,
                                  n_outside = self.n_outside,
                                  reflection_threshold=0.5,
                                  transmission_threshold=0.5)
        return m
    
    def _faces_default(self):
        fl = super(SolidRetroreflector, self)._faces_default()
        m = self.material
        fl.faces.append(CircularFace(owner=self, material=m, 
                                     z_plane=self.thickness,
                                     invert_normal=True))
        print("adding circular face", m)
        return fl
    
    def _thickness_changed(self, new_t):
        face = self.faces.faces[-1]
        print("setting thickness", face, face.z_plane)
        face.z_plane = new_t
        thickness = new_t
        size = max(thickness, self.diameter)*2
        cube = self.cube
        cube.set_bounds(0,size,0,size,0,size)
        
        cyl=self.cyl
        cyl.height = thickness
        cyl.center = (0, thickness/2,0)
        self.update=True
        
    def _diameter_changed(self, d):
        self.cyl.radius = d/2.
        size = max(self.thickness, d)*2
        cube = self.cube
        cube.set_bounds(0,size,0,size,0,size)
        self.update=True
        
    def _pipeline_default(self):
        cyl = self.cyl
        cyl.resolution = 63 #use a prime no to avoid vertex degeneracy
        # Important to ensure the cylinder end doesn't coincide with the corner of the cube
        cyl.height = self.thickness-1
        cyl.radius = self.diameter/2.0
        cyl.center = (0, (self.thickness/2)-1,0)
        print("thickness", self.thickness)
        print("diameter", self.diameter)
        size = max(self.thickness, self.diameter)*2
        cube = self.cube
        cube.set_bounds(0,size,0,size,0,size)
        
        tf = tvtk.Transform()
        tf.post_multiply()
        tf.rotate_x(-45.0)
        tf.rotate_wxyz(35.26438968275, 0,0,1)
        
        tfilt = tvtk.TransformFilter(input_connection=cube.output_port)
        tfilt.transform=tf
        
        tri1 = tvtk.TriangleFilter(input_connection=cyl.output_port)
        tri2 = tvtk.TriangleFilter(input_connection=tfilt.output_port)
        tri1.update()
        tri2.update()
        
        intersect = tvtk.BooleanOperationPolyDataFilter()
        intersect.operation = "intersection"
        intersect.add_input_connection(0, tri1.output_port)
        intersect.add_input_connection(1, tri2.output_port)
        intersect.tolerance = 1e-8
        
        tf2=tvtk.Transform()
        tf2.rotate_x(90.0)
        tf2.rotate_y(60.0)
        orient = tvtk.TransformFilter(input_connection=intersect.output_port,
                                      transform=tf2)
        
        norm = tvtk.PolyDataNormals(input_connection=orient.output_port)
        transF = tvtk.TransformFilter(input_connection=norm.output_port, 
                                      transform=self.transform)
        self.config_pipeline()
        return transF
    
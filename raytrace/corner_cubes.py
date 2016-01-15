#    Copyright 2009, Teraview Ltd., Bryan Cole
#
#    This file is part of Raytrace.
#
#    Raytrace is free software: you can redistribute it and/or modify
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
from itertools import chain, izip, islice, tee

#from raytrace.tracer import Optic, VTKOptic, normaliseVector, RaySegment,\
#             Traceable, NumEditor, dotprod, transformPoints, transformNormals

from raytrace.bases import Optic, Traceable, NumEditor, RaytraceObject
from raytrace.cfaces import ElipticalPlaneFace, CircularFace
from raytrace.ctracer import FaceList
from raytrace.mirrors import BaseMirror
from raytrace.cmaterials import PECMaterial


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
        append.add_input(p1.output)
        append.add_input(p2.output)
        append.add_input(p3.output)
        
        scale_f = tvtk.TransformFilter(input=append.output, transform=self.scale)
        
        trans = tvtk.Transform()
        trans.rotate_x(90.0)
        self.imp_cyl.transform = trans
        
        clip = tvtk.ClipPolyData(input=scale_f.output, clip_function = self.imp_cyl,
                                 inside_out=True)
                               
        transF2 = tvtk.TransformFilter(input=clip.output, transform=self.transform)
        self.config_pipeline()
        return transF2
    
    
class SolidRetroreflector(Optic, HollowRetroreflector):
    abstract = False
    yaml_tag = u"!SolidRetroreflector"
    n_inside = 1.5
    name = "Solid Retroreflector"
    
    offset = 0.0
    
    thickness = Float(20.0, desc="distance from input face to apex")
    
    vtk_grid = Instance(tvtk.ProgrammableSource, (), transient=True)
    cyl = Instance(tvtk.Cylinder, (), transient=True)
    
    traits_view = View(VGroup(
                       Traceable.uigroup,
                       Item('diameter', editor=NumEditor),
                       Item('thickness', editor=NumEditor),
                       Item('n_inside'),
                        ),
                   )
    
    
    def _faces_default(self):
        fl = super(SolidRetroreflector, self)._faces_default()
        m = self.material
        fl.faces.append(CircularFace(owner=self, material=m, 
                                     z_plane=self.thickness,
                                     invert_normal=True))
        print "adding circular face", m
        return fl
    
    def _thickness_changed(self, new_t):
        face = self.faces.faces[-1]
        print "setting thickness", face, face.z_plane
        face.z_plane = new_t
        self.vtk_grid.modified()
        self.update=True
        
    def _diameter_changed(self, d):
        self.vtk_grid.modified()
        self.cyl.radius = d/2.
        self.update=True
        
    def create_grid(self):
        r = 1.1*self.diameter/2.
        xmin, xmax = -r,r
        ymin, ymax = -r,r
        zmin, zmax = 0, self.thickness
        
        source = self.vtk_grid
        sp = source.structured_points_output
        size = 50
        dx = (xmax-xmin) / (size-1)
        dy = (ymax-ymin) / (size-1)
        dz = (zmax-zmin) / (size-1)
        sp.dimensions = (size,size,size)
        sp.whole_extent=(0,size,0,size,0,size)
        sp.origin = (xmin, ymin, zmin)
        sp.spacing = (dx, dy, dz)
        sp.set_update_extent_to_whole_extent()
        
    def _pipeline_default(self):
        grid = self.vtk_grid
        grid.set_execute_method(self.create_grid)
        grid.structured_points_output
        grid.modified()
        
        
        sqrt2 = numpy.sqrt(2)
        x = sqrt2*numpy.cos(numpy.pi*2./3)
        y = sqrt2*numpy.sin(numpy.pi*2./3)
        
        t = self.thickness
        planes = tvtk.Planes(points=[[0,0,0],[0,0,0],[0,0,0]],
                             normals=[[sqrt2,0,-1],[x,y,-1],[x,-y,-1]])
        tr = tvtk.Transform()
        tr.rotate_x(90.0)
        cyl = self.cyl
        cyl.transform = tr
        cyl.radius = self.diameter/2.
        union = tvtk.ImplicitBoolean(operation_type="intersection")
        union.add_function(planes)
        union.add_function(cyl)
        
        clip = tvtk.ClipVolume(input_connection=grid.output_port,
                                 clip_function=union,
                                 inside_out=True,
                                 mixed3d_cell_generation=True)
        
        topoly = tvtk.GeometryFilter(input_connection=clip.output_port)
        
        norm = tvtk.PolyDataNormals(input_connection=topoly.output_port)
        transF = tvtk.TransformFilter(input_connection=norm.output_port, 
                                      transform=self.transform)
        self.config_pipeline()
        grid.modified()
        return transF
    
    def _material_default(self):
        return Optic._material_default(self)